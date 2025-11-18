// ont_single.cpp - ONT single-read mapping pipeline
#include "ont_mapper.h"
#include "../../core/dna_utils.h"
#include <algorithm>

namespace rnamapper {

std::vector<Alignment> ONTMapper::map_read_single(
    const std::string &qname,
    const std::string &seq)
{
    ST.total_reads++;

    std::vector<Alignment> hits;
    const int L = (int)seq.size();

    // Check minimum length
    if (L < (int)IX.k) {
        ST.f_len_lt_k++;
        return hits;
    }

    // Prepare sequences
    std::string fwd = upper(seq);
    std::string rev = revcomp(fwd);
    auto f2b = pack2(fwd);
    auto r2b = pack2(rev);

    // Discover forward strand
    auto F = discover_ont(fwd, false);
    if (F.all_seeds_filtered) {
        ST.f_all_seeds_filtered++;
    }
    if (!F.shortlist.empty()) {
        ST.reads_with_shortlist++;
    }

    // Verify forward candidates
    int budgetF = choose_verify_budget_ont(F.shortlist, F.best_votes, F.best_runner);
    for (int i = 0; i < budgetF; i++) {
        auto a = verify_one_ont(fwd, f2b, F.shortlist[i], false);
        if (a) {
            hits.push_back(*a);
        }
    }

    // Discover reverse strand
    ST.reverse_attempted++;
    auto R = discover_ont(rev, true);
    if (R.all_seeds_filtered) {
        ST.f_all_seeds_filtered++;
    }

    // Verify reverse candidates
    int budgetR = choose_verify_budget_ont(R.shortlist, R.best_votes, R.best_runner);
    for (int i = 0; i < budgetR; i++) {
        auto a = verify_one_ont(rev, r2b, R.shortlist[i], true);
        if (a) {
            hits.push_back(*a);
        }
    }

    if (!F.shortlist.empty() || !R.shortlist.empty()) {
        ST.reads_with_shortlist++;
    }

    // Sort and deduplicate
    std::sort(hits.begin(), hits.end(), [](const Alignment &a, const Alignment &b) {
        if (a.score != b.score) return a.score > b.score;
        return a.mm < b.mm;
    });

    hits.erase(
        std::unique(hits.begin(), hits.end(), [](const Alignment &a, const Alignment &b) {
            return a.tid == b.tid && a.pos == b.pos && a.is_reverse == b.is_reverse;
        }),
        hits.end()
    );

    // Update stats
    if (hits.empty()) {
        if (F.shortlist.empty()) {
            ST.f_discover_no_shortlist++;
        } else {
            ST.f_verify_failed++;
        }
    } else {
        ST.reads_with_verified_hit++;
        ST.ok_mapped_reads++;
    }

    return hits;
}

} // namespace rnamapper
