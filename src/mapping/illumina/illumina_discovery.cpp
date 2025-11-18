// illumina_discovery.cpp - Discovery phase implementation
#include "illumina_mapper.h"
#include <algorithm>

namespace rnamapper {

IlluminaMapper::DiscoverRes IlluminaMapper::discover(const std::string &seq, bool rev) {
    next_read_epoch();
    DiscoverRes R;
    auto seeds = make_seeds(seq);
    R.had_any_seed = !seeds.empty();

    uint32_t posts_seen = 0;
    uint16_t best_votes = 0, best_runner = 0;
    uint32_t best_tid = UINT32_MAX;
    const uint32_t KEEP_TOP_TIDS = 32;
    tid_list.clear();

    auto mark_tid = [&](uint32_t tid, uint32_t add) {
        if (sieve_gen[tid] != per_tid_gen) {
            sieve_gen[tid] = per_tid_gen;
            sieve_cnt[tid] = (uint16_t)std::min<uint32_t>(65535u, add);
            tid_list.push_back(tid);
        } else {
            uint32_t v = (uint32_t)sieve_cnt[tid] + add;
            sieve_cnt[tid] = (uint16_t)std::min<uint32_t>(65535u, v);
        }
    };

    bool saw_effective_posting = false;
    int processed = 0;

    int s1;
    if (sieve_seeds <= 0)
        s1 = 0;
    else if ((int)seeds.size() <= 1)
        s1 = 0;
    else
        s1 = std::min<int>(sieve_seeds, (int)seeds.size() - 1);

    for (int i = 0; i < s1; i++) {
        const auto &s = seeds[i];
        if (skip_common && s.common) continue;
        if (s.df_bytes >= df1_bytes) continue;
        IX.for_each_posting_breakable(s.key_idx, [&](uint32_t tid, uint32_t) {
            posts_seen++;
            mark_tid(tid, 1);
            saw_effective_posting = true;
            return true;
        }, prefetch);
    }

    std::vector<uint8_t> allow;
    bool have_allow = false;
    if (s1 > 0 && !tid_list.empty()) {
        struct Pair { uint16_t c; uint32_t tid; };
        auto cmp_min = [](const Pair &a, const Pair &b) { return a.c > b.c; };
        std::vector<Pair> heap;
        heap.reserve(KEEP_TOP_TIDS + 8);
        for (uint32_t tid : tid_list) {
            uint16_t c = sieve_cnt[tid];
            if ((int)heap.size() < (int)KEEP_TOP_TIDS) {
                heap.push_back({c, tid});
                std::push_heap(heap.begin(), heap.end(), cmp_min);
            } else if (c > heap.front().c) {
                std::pop_heap(heap.begin(), heap.end(), cmp_min);
                heap.back() = {c, tid};
                std::push_heap(heap.begin(), heap.end(), cmp_min);
            }
        }
        std::sort(heap.begin(), heap.end(), [](const Pair &a, const Pair &b) {
            return a.c > b.c;
        });
        allow.assign((IX.n_targets + 7) >> 3, 0);
        auto set_allow = [&](uint32_t tid) {
            allow[tid >> 3] |= (uint8_t)(1u << (tid & 7));
        };
        for (const auto &p : heap) set_allow(p.tid);
        have_allow = !heap.empty();
    }

    for (size_t i = s1; i < seeds.size(); ++i) {
        const auto &s = seeds[i];
        if (skip_common && s.common) continue;
        if (s.df_bytes >= df2_bytes) continue;

        IX.for_each_posting_breakable(s.key_idx, [&](uint32_t tid, uint32_t pos) -> bool {
            posts_seen++;
            if (have_allow) {
                if (((allow[tid >> 3] >> (tid & 7)) & 1u) == 0) {
                    if (post_budget && posts_seen >= post_budget && best_tid != UINT32_MAX) return false;
                    return true;
                }
            }
            saw_effective_posting = true;

            int32_t off = (int32_t)pos - (int32_t)s.pos;
            int32_t off_q = off >> bin_shift;
            uint8_t flags = 0;
            if (IX.is_jx_tid(tid)) {
                uint32_t jxid = IX.jx_index(tid);
                const auto &J = IX.jx[jxid];
                uint32_t mid = J.len_bp >> 1;
                flags |= (pos < mid) ? 1 : 2;
            }
            uint16_t c = bump_hh(tid, off_q, flags);
            on_vote(tid, off_q, c);

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

            if ((processed >= min_seeds_before_early) &&
                best_tid != UINT32_MAX &&
                ((int)best_votes - (int)best_runner) >= lead_margin) {
                return false;
            }
            if (post_budget && posts_seen >= post_budget && best_tid != UINT32_MAX) return false;
            return true;
        }, prefetch);

        processed++;
        int need = min_votes;
        if (best_tid != UINT32_MAX &&
            best_votes >= need &&
            (int)best_votes - (int)best_runner >= lead_margin &&
            processed >= min_seeds_before_early) {
            break;
        }
        if (post_budget && posts_seen >= post_budget && best_tid != UINT32_MAX) break;
    }

    std::vector<Cand> tmp;
    tmp.reserve(hh_touched.size());
    for (uint32_t tid : hh_touched) {
        auto &tb = per_tid[tid];
        if (tb.gen != per_tid_gen) continue;
        const auto &H = hh[tid];
        int win = -1;
        for (int i = 0; i < 4; i++) {
            if (H.cnt[i] == 0) continue;
            if (H.off_q[i] == tb.best_off_q && H.cnt[i] == tb.best) {
                win = i;
                break;
            }
        }
        if (win == -1) {
            for (int i = 0; i < 4; i++)
                if (H.cnt[i] == tb.best) {
                    win = i;
                    break;
                }
        }
        if (win != -1) {
            tmp.push_back({tid, H.off_q[win], H.cnt[win], H.flg[win]});
        }
    }

    tmp.erase(std::remove_if(tmp.begin(), tmp.end(), [&](const Cand &c) {
        if (!IX.is_jx_tid(c.tid)) return false;
        return ((c.flags & 3) != 3);
    }), tmp.end());

    std::sort(tmp.begin(), tmp.end(), [](const Cand &a, const Cand &b) {
        if (a.votes != b.votes) return a.votes > b.votes;
        if (a.tid != b.tid) return a.tid < b.tid;
        return a.off_q < b.off_q;
    });

    R.shortlist = std::move(tmp);
    R.best_votes = best_votes;
    R.best_runner = best_runner;

    if (R.had_any_seed && !saw_effective_posting)
        R.all_seeds_filtered = true;

    return R;
}

} // namespace rnamapper
