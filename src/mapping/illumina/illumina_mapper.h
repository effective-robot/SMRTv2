// illumina_mapper.h - Illumina paired-end RNA-seq mapper
#pragma once

#include "../mapper_base.h"
#include "../../core/types.h"
#include "../../core/dna_utils.h"
#include "../../index/index.h"
#include "../../alignment/hamming.h"
#include <vector>
#include <optional>
#include <unordered_set>
#include <cstdint>
#include <climits>

namespace rnamapper {

class IlluminaMapper : public MapperBase {
public:

    // Algorithm components
    Hamming2B ham;

    // Tunable parameters
    int stride = 1, bin_shift = 3, min_votes = 3, lead_margin = 2, min_seeds_before_early = 6;
    bool skip_common = true;
    int verify_topk = 12, verify_topk_lo = 5;
    uint32_t df1_bytes = 128, df2_bytes = 512, post_budget = 0;
    bool prefetch = true, dedup_seeds = false;
    int sieve_seeds = 10;
    bool skip_rev_on_clear = true;
    int min_ins = 100, max_ins = 800;
    bool output_nm_md = false;

    // Internal data structures
    struct HH4 {
        int32_t off_q[4];
        uint16_t cnt[4];
        uint8_t flg[4];
        uint32_t gen;
        HH4();
    };

    std::vector<HH4> hh;
    std::vector<uint32_t> hh_touched;

    struct TBest {
        uint16_t best = 0, second = 0;
        uint32_t gen = 0;
        int32_t best_off_q = INT32_MIN;
    };

    std::vector<TBest> per_tid;
    uint32_t per_tid_gen = 1;
    std::vector<uint16_t> sieve_cnt;
    std::vector<uint32_t> sieve_gen;
    std::vector<uint32_t> tid_list;

    // Seed structure
    struct Seed {
        uint64_t kmer;
        int pos;
        uint32_t key_idx, df_bytes;
        bool common;
    };

    // Candidate structure
    struct Cand {
        uint32_t tid;
        int32_t off_q;
        uint16_t votes;
        uint8_t flags;
    };

    // Discovery result
    struct DiscoverRes {
        std::vector<Cand> shortlist;
        uint16_t best_votes = 0, best_runner = 0;
        bool had_any_seed = false;
        bool all_seeds_filtered = false;
    };

    explicit IlluminaMapper(const IndexVX &ix, FailureStats &st);

    // Core mapping functions
    void next_read_epoch() override;
    std::vector<Seed> make_seeds(const std::string &seq) const;
    DiscoverRes discover(const std::string &seq, bool rev);
    std::optional<Alignment> verify_one(const std::string &read,
                                       const std::vector<uint8_t> &read2b,
                                       const Cand &c,
                                       bool strand_is_rev,
                                       double mm_frac_tx = 0.10,
                                       double mm_frac_jx = 0.05) const;

    std::vector<Alignment> map_read_single(const std::string &qname, const std::string &seq) override;

    // Paired-end mapping
    void map_paired(const std::string &r1, const std::string &r2,
                   const std::string &out_path, uint64_t max_pairs = 0, bool do_rescue = false) override;

    void map_single_end(const std::string &reads, const std::string &out_path,
                       uint64_t max_reads = 0, bool do_rescue = false) override;

    // Helper functions (declared inline for performance but implemented in .cpp)
    uint16_t bump_hh(uint32_t tid, int32_t off_q, uint8_t flags);
    void on_vote(uint32_t tid, int32_t off_q, uint16_t cnt);
    int choose_verify_budget(const std::vector<Cand> &C, uint16_t best_votes, uint16_t best_runner) const;
    bool forward_is_clear(uint16_t best_votes, uint16_t best_runner) const;
    bool tid_compatible_same_gene(uint32_t a, uint32_t b) const;

    // Genomic projection
    GenomicPos transcript_to_genomic(uint32_t tid, int tx_pos, int read_len) const override;

    // MAPQ calculation
    uint8_t calculate_mapq(const Alignment *best, const std::vector<Alignment> &all_hits, int read_len) const override;

    // Insert size accessors
    int get_min_ins() const override { return min_ins; }
    int get_max_ins() const override { return max_ins; }

    // Pair scoring
    struct PairScore {
        int r1_score;
        int r2_score;
        int concordance_bonus;
        int insert_quality;
        int total;
        bool operator<(const PairScore &other) const { return total > other.total; }
    };

    PairScore score_pair(const Alignment &a1, const Alignment &a2,
                        const ReadRecord &rec1, const ReadRecord &rec2) const;

    struct BestPair {
        Alignment *best1;
        Alignment *best2;
        PairScore score;
        bool found;
    };

    BestPair select_best_pair(std::vector<Alignment> &H1, std::vector<Alignment> &H2,
                             const ReadRecord &rec1, const ReadRecord &rec2) const;

    bool attempt_mate_rescue(std::vector<Alignment> &H1, Alignment *&best1,
                            std::vector<Alignment> &H2, Alignment *&best2,
                            const ReadRecord &rec1, const ReadRecord &rec2);

    // Rescue/reprocess
    struct ReprocessResult {
        bool recovered = false;
        Alignment aln;
        std::string tier_name;
        std::string notes;
    };

    ReprocessResult attempt_reprocess_single(const std::string &qname, const std::string &seq,
                                            const Alignment *orig_aln = nullptr);

private:
    static std::string strip_isoform(const std::string &s);
};

} // namespace rnamapper
