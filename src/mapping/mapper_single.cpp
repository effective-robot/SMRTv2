// mapper_single.cpp - Single-end read mapping
#include "mapper.h"
#include <algorithm>

namespace rnamapper {

std::vector<Alignment> Mapper::map_read_single(const std::string &qname, const std::string &seq) {
    ST.total_reads++;
    std::vector<Alignment> hits;
    const int L = (int)seq.size();
    if (L < (int)IX.k) {
        ST.f_len_lt_k++;
        return hits;
    }

    std::string fwd = upper(seq);
    std::string rev = revcomp(fwd);
    auto f2b = pack2(fwd);
    auto r2b = pack2(rev);

    auto F = discover(fwd, false);
    if (F.all_seeds_filtered) ST.f_all_seeds_filtered++;
    if (!F.shortlist.empty()) ST.reads_with_shortlist++;

    int budgetF = choose_verify_budget(F.shortlist, F.best_votes, F.best_runner);
    for (int i = 0; i < budgetF; i++) {
        auto a = verify_one(fwd, f2b, F.shortlist[i], false);
        if (a) hits.push_back(*a);
    }

    bool do_reverse = true;
    if (skip_rev_on_clear && forward_is_clear(F.best_votes, F.best_runner) && !hits.empty()) {
        do_reverse = false;
        ST.reverse_skipped_clear++;
    }

    if (do_reverse) {
        ST.reverse_attempted++;
        auto R = discover(rev, true);
        if (R.all_seeds_filtered) ST.f_all_seeds_filtered++;

        int budgetR = choose_verify_budget(R.shortlist, R.best_votes, R.best_runner);
        for (int i = 0; i < budgetR; i++) {
            auto a = verify_one(rev, r2b, R.shortlist[i], true);
            if (a) hits.push_back(*a);
        }

        if (!F.shortlist.empty() || !R.shortlist.empty())
            ST.reads_with_shortlist++;
    }

    std::sort(hits.begin(), hits.end(), [](const Alignment &a, const Alignment &b) {
        if (a.score != b.score) return a.score > b.score;
        return a.mm < b.mm;
    });
    hits.erase(std::unique(hits.begin(), hits.end(), [](const Alignment &a, const Alignment &b) {
        return a.tid == b.tid && a.pos == b.pos && a.is_reverse == b.is_reverse;
    }), hits.end());

    if (hits.empty()) {
        bool f_empty = F.shortlist.empty();
        if (f_empty) ST.f_discover_no_shortlist++;
        else ST.f_verify_failed++;
    } else {
        ST.reads_with_verified_hit++;
    }

    // Tie detection logic
    if (hits.size() >= 2) {
        int score_diff = hits[0].score - hits[1].score;
        GenomicPos gp1 = transcript_to_genomic(hits[0].tid, hits[0].pos, (int)seq.size());
        GenomicPos gp2 = transcript_to_genomic(hits[1].tid, hits[1].pos, (int)seq.size());

        if (score_diff <= 2 && gp1.valid && gp2.valid && gp1.chrom != gp2.chrom) {
            hits[0].is_ambiguous = true;
        }
    } else if (hits.size() == 1) {
        hits[0].is_ambiguous = true;
    }

    return hits;
}

} // namespace rnamapper
