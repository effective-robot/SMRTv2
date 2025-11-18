// ont_discovery.cpp - ONT candidate discovery implementation
#include "ont_mapper.h"
#include <algorithm>

namespace rnamapper {

void ONTMapper::next_read_epoch() {
    per_tid_gen++;
    if (per_tid_gen == 0) {
        per_tid_gen = 1;
        // Reset all structures
        for (auto &t : per_tid) {
            t.gen = 0;
            t.best = t.second = 0;
            t.best_off_q = INT32_MIN;
        }
        for (auto &h : hh) {
            h.gen = 0;
            for (int i = 0; i < 4; i++) {
                h.off_q[i] = INT32_MIN;
                h.cnt[i] = 0;
                h.flg[i] = 0;
            }
        }
    }
    hh_touched.clear();
}

uint16_t ONTMapper::bump_hh(uint32_t tid, int32_t off_q, uint8_t flags) {
    HH4 &H = hh[tid];
    if (H.gen != per_tid_gen) {
        H.gen = per_tid_gen;
        for (int i = 0; i < 4; ++i) {
            H.off_q[i] = INT32_MIN;
            H.cnt[i] = 0;
            H.flg[i] = 0;
        }
    }

    // Try to find existing offset
    for (int i = 0; i < 4; ++i) {
        if (H.off_q[i] == off_q) {
            H.cnt[i] += 1;
            H.flg[i] |= flags;
            return H.cnt[i];
        }
    }

    // Find empty slot
    for (int i = 0; i < 4; ++i) {
        if (H.off_q[i] == INT32_MIN) {
            H.off_q[i] = off_q;
            H.cnt[i] = 1;
            H.flg[i] = flags;
            hh_touched.push_back(tid);
            return 1;
        }
    }

    // Replace worst
    int worst = 0;
    for (int i = 1; i < 4; ++i) {
        if (H.cnt[i] < H.cnt[worst]) worst = i;
    }
    H.off_q[worst] = off_q;
    H.cnt[worst] = 1;
    H.flg[worst] = flags;
    hh_touched.push_back(tid);
    return 1;
}

void ONTMapper::on_vote(uint32_t tid, int32_t off_q, uint16_t cnt) {
    auto &t = per_tid[tid];
    if (t.gen != per_tid_gen) {
        t.gen = per_tid_gen;
        t.best = cnt;
        t.second = 0;
        t.best_off_q = off_q;
        return;
    }

    if (cnt >= t.best) {
        if (off_q != t.best_off_q) {
            t.second = t.best;
        }
        t.best = cnt;
        t.best_off_q = off_q;
    } else if (cnt > t.second) {
        t.second = cnt;
    }
}

ONTMapper::DiscoverRes ONTMapper::discover_ont(const std::string &seq, bool is_rev) {
    next_read_epoch();
    DiscoverRes R;

    // Generate seeds
    auto seeds = make_seeds_ont(seq);
    R.had_any_seed = !seeds.empty();

    if (seeds.empty()) {
        return R;
    }

    // Process seeds and aggregate votes
    uint16_t best_votes = 0, best_runner = 0;
    uint32_t best_tid = UINT32_MAX;

    // Use relaxed DF filtering for ONT
    uint32_t df_limit = 2048;  // Higher than Illumina's 512

    bool saw_effective_posting = false;

    for (const auto &s : seeds) {
        // Skip high-DF seeds
        if (s.df_bytes >= df_limit) {
            continue;
        }

        // Process postings for this seed
        IX.for_each_posting_breakable(
            s.key_idx,
            [&](uint32_t tid, uint32_t pos) -> bool {
                saw_effective_posting = true;

                // Compute quantized offset
                int32_t off = (int32_t)pos - (int32_t)s.pos;
                int32_t off_q = off >> bin_shift;

                // Junction flags
                uint8_t flags = 0;
                if (IX.is_jx_tid(tid)) {
                    uint32_t jxid = IX.jx_index(tid);
                    const auto &J = IX.jx[jxid];
                    uint32_t mid = J.len_bp >> 1;
                    flags |= (pos < mid) ? 1 : 2;
                }

                // Update votes
                uint16_t c = bump_hh(tid, off_q, flags);
                on_vote(tid, off_q, c);

                // Track best
                const auto &tb = per_tid[tid];
                if (tb.gen == per_tid_gen && tb.best == c) {
                    if (c > best_votes) {
                        best_votes = c;
                        best_tid = tid;
                        best_runner = tb.second;
                    } else if (tid == best_tid) {
                        best_runner = tb.second;
                    }
                }

                return true;  // Continue processing
            },
            true  // prefetch
        );
    }

    // Build shortlist from touched transcripts
    std::vector<Cand> tmp;
    tmp.reserve(hh_touched.size());

    for (uint32_t tid : hh_touched) {
        auto &tb = per_tid[tid];
        if (tb.gen != per_tid_gen) continue;

        const auto &H = hh[tid];
        int win = -1;

        // Find winning offset
        for (int i = 0; i < 4; i++) {
            if (H.cnt[i] == 0) continue;
            if (H.off_q[i] == tb.best_off_q && H.cnt[i] == tb.best) {
                win = i;
                break;
            }
        }

        if (win == -1) {
            for (int i = 0; i < 4; i++) {
                if (H.cnt[i] == tb.best) {
                    win = i;
                    break;
                }
            }
        }

        if (win != -1) {
            tmp.push_back({tid, H.off_q[win], H.cnt[win], H.flg[win]});
        }
    }

    // Filter junction candidates
    tmp.erase(
        std::remove_if(tmp.begin(), tmp.end(), [&](const Cand &c) {
            if (!IX.is_jx_tid(c.tid)) return false;
            return ((c.flags & 3) != 3);  // Need both sides
        }),
        tmp.end()
    );

    // Sort by votes
    std::sort(tmp.begin(), tmp.end(), [](const Cand &a, const Cand &b) {
        if (a.votes != b.votes) return a.votes > b.votes;
        if (a.tid != b.tid) return a.tid < b.tid;
        return a.off_q < b.off_q;
    });

    R.shortlist = std::move(tmp);
    R.best_votes = best_votes;
    R.best_runner = best_runner;

    if (R.had_any_seed && !saw_effective_posting) {
        R.all_seeds_filtered = true;
    }

    return R;
}

} // namespace rnamapper
