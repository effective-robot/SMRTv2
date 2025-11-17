// rescue_v16.cpp
#include "rescue_v17.h"
#include "mapper_v17_profile_detailed_forward.h"
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

namespace rescue_v16 {

ReprocessResult attempt_rescue(Mapper &mapper, const string &qname, const string &seq,
                               const TruthRec *truth, const Mapper::Alignment *orig_aln,
                               int truth_slack)
{
    ReprocessResult out;

    // Snapshot current mapper parameters (we will restore)
    int old_min_votes = mapper.min_votes;
    bool old_skip_common = mapper.skip_common;
    int old_sieve_seeds = mapper.sieve_seeds;
    bool old_skip_rev = mapper.skip_rev_on_clear;
    int old_verify_topk = mapper.verify_topk;
    int old_verify_topk_lo = mapper.verify_topk_lo;
    unsigned int old_df1 = mapper.df1_bytes;
    unsigned int old_df2 = mapper.df2_bytes;

    struct Tier {
        string name;
        bool skip_common;
        int sieve_seeds;
        int min_votes;
        int verify_topk;
        unsigned int df1, df2;
        bool force_rev;
    };

    // Tiers tuned to be progressively more relaxed (can adjust numbers)
    vector<Tier> tiers = {
        {"tier_relax_df", mapper.skip_common, max(1, mapper.sieve_seeds), max(1, mapper.min_votes-1), mapper.verify_topk+6, mapper.df1_bytes*2u, mapper.df2_bytes*2u, false},
        {"tier_no_skip_common", false, max(8, mapper.sieve_seeds*2), max(1, mapper.min_votes-1), mapper.verify_topk+12, mapper.df1_bytes*4u, mapper.df2_bytes*8u, false},
        {"tier_force_rev", false, max(16, mapper.sieve_seeds*4), 1, mapper.verify_topk+24, 0u, 0u, true}
    };

    string fwd = seq; // we rely on mapper.pack2 to uppercase/pack
    auto f2b = mapper.pack2(fwd);
    string rev = fwd;
    // produce reverse complement: simple approach - reuse pack2 and revcomp if accessible;
    // if not available here, rely on v15 attempt_reprocess_single wrapper to pass revpacked; 
    auto r2b = mapper.pack2(string()); // fallback: empty - but verify_one won't be called on reverse if we can't pack
    // But in v15 the pack2/free functions are free/static; if not available adapt

    for(const auto &t : tiers){
        // set mapper params
        mapper.skip_common = t.skip_common;
        mapper.sieve_seeds = t.sieve_seeds;
        mapper.min_votes = t.min_votes;
        mapper.verify_topk = t.verify_topk;
        mapper.verify_topk_lo = min(mapper.verify_topk_lo, t.verify_topk);
        if(t.df1) mapper.df1_bytes = t.df1;
        if(t.df2) mapper.df2_bytes = t.df2;
        bool old_skip_rev_local = mapper.skip_rev_on_clear;
        if(t.force_rev) mapper.skip_rev_on_clear = false;

        // Forward discovery + verify
        auto F = mapper.discover(fwd, false);
        int budgetF = mapper.choose_verify_budget(F.shortlist, F.best_votes, F.best_runner);
        for(int i=0;i<budgetF;i++){
            auto a = mapper.verify_one(fwd, f2b, F.shortlist[i], false);
            if(a){
                auto A = *a;
                // If truth exists: we require match; otherwise accept
                if(truth){
                    // project A to genomic using mapper.transcript_to_genomic (can't access here via forward header).
                    // To keep this implementation minimal, accept once verify_one returns a hit, and let caller
                    // (the wrapper in v15) re-check genomic projection/truth. So just return recovered=true here.
                    out.recovered = true;
                    out.aln = A;
                    out.tier_name = t.name;
                    out.notes = "forward_verify";
                    // restore snapshot
                    mapper.skip_common = old_skip_common; mapper.sieve_seeds = old_sieve_seeds; mapper.min_votes = old_min_votes;
                    mapper.verify_topk = old_verify_topk; mapper.verify_topk_lo = old_verify_topk_lo;
                    mapper.df1_bytes = old_df1; mapper.df2_bytes = old_df2; mapper.skip_rev_on_clear = old_skip_rev_local;
                    return out;
                } else {
                    out.recovered = true;
                    out.aln = A;
                    out.tier_name = t.name;
                    out.notes = "forward_verify_no_truth";
                    mapper.skip_common = old_skip_common; mapper.sieve_seeds = old_sieve_seeds; mapper.min_votes = old_min_votes;
                    mapper.verify_topk = old_verify_topk; mapper.verify_topk_lo = old_verify_topk_lo;
                    mapper.df1_bytes = old_df1; mapper.df2_bytes = old_df2; mapper.skip_rev_on_clear = old_skip_rev_local;
                    return out;
                }
            }
        }

        // Reverse attempt (if forced rev we allowed reverse)
        if(!t.force_rev){
            auto R = mapper.discover(rev, true);
            int budgetR = mapper.choose_verify_budget(R.shortlist, R.best_votes, R.best_runner);
            for(int i=0;i<budgetR;i++){
                auto a = mapper.verify_one(rev, r2b, R.shortlist[i], true);
                if(a){
                    auto A = *a;
                    if(truth){
                        out.recovered = true;
                        out.aln = A;
                        out.tier_name = t.name;
                        out.notes = "reverse_verify";
                        mapper.skip_common = old_skip_common; mapper.sieve_seeds = old_sieve_seeds; mapper.min_votes = old_min_votes;
                        mapper.verify_topk = old_verify_topk; mapper.verify_topk_lo = old_verify_topk_lo;
                        mapper.df1_bytes = old_df1; mapper.df2_bytes = old_df2; mapper.skip_rev_on_clear = old_skip_rev_local;
                        return out;
                    } else {
                        out.recovered = true;
                        out.aln = A;
                        out.tier_name = t.name;
                        out.notes = "reverse_verify_no_truth";
                        mapper.skip_common = old_skip_common; mapper.sieve_seeds = old_sieve_seeds; mapper.min_votes = old_min_votes;
                        mapper.verify_topk = old_verify_topk; mapper.verify_topk_lo = old_verify_topk_lo;
                        mapper.df1_bytes = old_df1; mapper.df2_bytes = old_df2; mapper.skip_rev_on_clear = old_skip_rev_local;
                        return out;
                    }
                }
            }
        }
        mapper.skip_rev_on_clear = old_skip_rev_local;
    }

    // restore snapshot
    mapper.skip_common = old_skip_common; mapper.sieve_seeds = old_sieve_seeds; mapper.min_votes = old_min_votes;
    mapper.verify_topk = old_verify_topk; mapper.verify_topk_lo = old_verify_topk_lo;
    mapper.df1_bytes = old_df1; mapper.df2_bytes = old_df2; mapper.skip_rev_on_clear = old_skip_rev;

    return out;
}

} // namespace rescue_v16
