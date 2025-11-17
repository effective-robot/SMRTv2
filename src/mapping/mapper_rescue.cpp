// mapper_rescue.cpp - Rescue/reprocess logic for unmapped reads
#include "mapper.h"
#include <unordered_set>
#include <sstream>
#include <algorithm>

namespace rnamapper {

Mapper::ReprocessResult Mapper::attempt_reprocess_single(const std::string &qname, const std::string &seq,
                                                         const Alignment *orig_aln) {
    ReprocessResult out;

    // Snapshot mapper state to restore later
    int old_min_votes = min_votes;
    int old_lead_margin = lead_margin;
    bool old_skip_common = skip_common;
    int old_sieve_seeds = sieve_seeds;
    bool old_skip_rev = skip_rev_on_clear;
    int old_verify_topk = verify_topk;
    int old_verify_topk_lo = verify_topk_lo;
    uint32_t old_df1 = df1_bytes;
    uint32_t old_df2 = df2_bytes;
    uint32_t old_post_budget = post_budget;
    int old_stride = stride;
    int old_bin_shift = bin_shift;
    bool old_dedup_seeds = dedup_seeds;

    ST.reprocess_attempts++;

    // Rescue tuning: dense seeding and relaxed controls
    const uint32_t DF_DISABLE = 1000000u;
    min_votes = 1;
    lead_margin = 0;
    skip_common = false;
    sieve_seeds = std::max(80, old_sieve_seeds * 4);
    df1_bytes = DF_DISABLE;
    df2_bytes = DF_DISABLE;
    post_budget = 0;
    skip_rev_on_clear = false;
    stride = 1;
    bin_shift = old_bin_shift;
    verify_topk = std::max(old_verify_topk, 200);
    verify_topk_lo = std::min(verify_topk, old_verify_topk_lo + 40);
    dedup_seeds = false;

    auto restore_state = [&]() {
        min_votes = old_min_votes;
        lead_margin = old_lead_margin;
        skip_common = old_skip_common;
        sieve_seeds = old_sieve_seeds;
        skip_rev_on_clear = old_skip_rev;
        verify_topk = old_verify_topk;
        verify_topk_lo = old_verify_topk_lo;
        df1_bytes = old_df1;
        df2_bytes = old_df2;
        post_budget = old_post_budget;
        stride = old_stride;
        bin_shift = old_bin_shift;
        dedup_seeds = old_dedup_seeds;
    };

    std::string fwd = upper(seq);
    std::string rev = revcomp(fwd);
    auto f2b = pack2(fwd);
    auto r2b = pack2(rev);

    // Try discover+verify helper
    auto try_discover_verify = [&](const std::string &s, const std::vector<uint8_t> &s2b, bool is_rev) -> std::optional<Alignment> {
        auto R = discover(s, is_rev);
        if (R.all_seeds_filtered) ST.f_all_seeds_filtered++;
        if (R.shortlist.empty()) return std::nullopt;

        int budget = choose_verify_budget(R.shortlist, R.best_votes, R.best_runner);
        bool has_clear = (R.best_votes >= (uint16_t)min_votes) && ((int)R.best_votes - (int)R.best_runner >= lead_margin);
        if (!has_clear) budget = std::min<int>((int)R.shortlist.size(), verify_topk);
        budget = std::max(1, std::min((int)R.shortlist.size(), budget));

        const double mm_frac_tx = 0.10;
        const double mm_frac_jx = 0.06;

        std::optional<Alignment> bestA = std::nullopt;
        for (int i = 0; i < budget; i++) {
            auto a = verify_one(s, s2b, R.shortlist[i], is_rev, mm_frac_tx, mm_frac_jx);
            if (!a) continue;
            Alignment A = *a;
            if (!bestA.has_value() || (A.score > bestA->score) || (A.score == bestA->score && A.mm < bestA->mm)) {
                bestA = A;
            }
            if (A.mm == 0) return bestA;
        }
        return bestA;
    };

    // Try forward
    if (auto Aopt = try_discover_verify(fwd, f2b, false); Aopt.has_value()) {
        Alignment A = *Aopt;
        out.recovered = true;
        out.aln = A;
        out.tier_name = "rescue_discover_forward";
        out.notes = "discover_verify";
        ST.reprocess_recovered++;
        restore_state();
        return out;
    }

    // Try reverse
    if (auto Aopt = try_discover_verify(rev, r2b, true); Aopt.has_value()) {
        Alignment A = *Aopt;
        out.recovered = true;
        out.aln = A;
        out.tier_name = "rescue_discover_reverse";
        out.notes = "discover_verify_rev";
        ST.reprocess_recovered++;
        restore_state();
        return out;
    }

    // Fallback: dense seed generation and manual aggregation
    struct SeedDiag {
        uint64_t kmer;
        int pos;
        bool found_key;
        uint32_t key_idx;
    };
    std::vector<SeedDiag> seeds;
    seeds.reserve(std::max(1, (int)fwd.size() - (int)IX.k + 1));
    const int K = (int)IX.k;
    uint64_t KMASK = (K == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * K)) - 1ULL);
    uint64_t roll = 0;
    int good = 0;

    for (size_t i = 0; i < fwd.size(); ++i) {
        int b = nt2b(fwd[i]);
        if (b < 0) {
            roll = 0;
            good = 0;
            continue;
        }
        roll = ((roll << 2) | (uint64_t)b) & KMASK;
        if (++good >= K) {
            int start = (int)i - K + 1;
            if ((start % std::max(1, stride)) != 0) continue;
            uint64_t kmer = roll;
            if (IX.uses_canonical) {
                uint64_t rc = 0, kx = kmer;
                for (int j = 0; j < K; j++) {
                    uint64_t base = kx & 3ULL;
                    rc = (rc << 2) | (3ULL - base);
                    kx >>= 2;
                }
                if (rc < kmer) kmer = rc;
            }
            auto [ok, idx] = IX.find_key(kmer);
            if (ok) seeds.push_back({kmer, start, true, idx});
        }
    }

    // Manual aggregation
    next_read_epoch();
    std::unordered_set<uint32_t> used_keyidxs;
    for (const auto &sd : seeds) {
        if (!sd.found_key) continue;
        uint32_t kid = sd.key_idx;
        if (!used_keyidxs.insert(kid).second) continue;
        IX.for_each_posting_breakable(kid, [&](uint32_t tid, uint32_t pos) -> bool {
            int32_t off_q = (int32_t)pos - (int32_t)sd.pos;
            off_q = off_q >> bin_shift;
            uint8_t flags = 0;
            if (IX.is_jx_tid(tid)) {
                uint32_t jxid = IX.jx_index(tid);
                const auto &J = IX.jx[jxid];
                uint32_t mid = J.len_bp >> 1;
                flags |= (pos < mid) ? 1 : 2;
            }
            uint16_t c = bump_hh(tid, off_q, flags);
            on_vote(tid, off_q, c);
            return true;
        }, prefetch);
    }

    // Collect shortlist
    std::vector<Cand> shortlist;
    shortlist.reserve(hh_touched.size());
    for (uint32_t tid : hh_touched) {
        auto &tb = per_tid[tid];
        if (tb.gen != per_tid_gen) continue;
        auto &H = hh[tid];
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
        if (win != -1) shortlist.push_back({tid, H.off_q[win], H.cnt[win], H.flg[win]});
    }

    shortlist.erase(std::remove_if(shortlist.begin(), shortlist.end(), [&](const Cand &c) {
        if (!IX.is_jx_tid(c.tid)) return false;
        return ((c.flags & 3) != 3);
    }), shortlist.end());

    std::sort(shortlist.begin(), shortlist.end(), [](const Cand &a, const Cand &b) {
        if (a.votes != b.votes) return a.votes > b.votes;
        if (a.tid != b.tid) return a.tid < b.tid;
        return a.off_q < b.off_q;
    });

    int budgetVerify = std::min<int>((int)shortlist.size(), verify_topk);
    const double mm_frac_tx = 0.10;
    const double mm_frac_jx = 0.06;
    std::optional<Alignment> bestA = std::nullopt;
    for (int i = 0; i < budgetVerify; i++) {
        auto a = verify_one(fwd, f2b, shortlist[i], false, mm_frac_tx, mm_frac_jx);
        if (!a) continue;
        Alignment A = *a;
        if (!bestA.has_value() || (A.score > bestA->score) || (A.score == bestA->score && A.mm < bestA->mm)) {
            bestA = A;
        }
        if (A.mm == 0) break;
    }

    if (bestA.has_value()) {
        Alignment A = *bestA;
        out.recovered = true;
        out.aln = A;
        out.tier_name = "rescue_postings_verify";
        out.notes = "fallback_verify";
        ST.reprocess_recovered++;
        restore_state();
        return out;
    }

    ST.reprocess_failed++;
    restore_state();
    return out;
}

} // namespace rnamapper
