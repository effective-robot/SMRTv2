// ont_mapper.h - ONT long-read RNA-seq mapper
#pragma once

#include "../mapper_base.h"
#include "../../core/types.h"
#include "../../index/index.h"
#include "../../alignment/hamming.h"  // Temporary - will replace with ONT aligner
#include <vector>
#include <optional>
#include <cstdint>

namespace rnamapper {

/**
 * ONTMapper - RNA-seq mapper optimized for ONT long reads
 *
 * Key differences from IlluminaMapper:
 * - Adaptive stride based on read length (50bp-22kb range)
 * - Sparse seeding for long reads (minimizers or sparse k-mers)
 * - Relaxed error tolerance (10-20% vs Illumina's <1%)
 * - Custom alignment for high error rate sequences
 * - No mate rescue (ONT is single-molecule)
 *
 * Implementation plan:
 * - Step 4: Adaptive seeding strategy
 * - Step 5: Candidate discovery with relaxed thresholds
 * - Step 6: ONT-specific verification
 * - Step 7: Single-read mapping pipeline
 * - Step 8: Genomic projection and MAPQ
 * - Step 9: File I/O and integration
 */
class ONTMapper : public MapperBase {
public:
    // Seed structure (similar to Illumina but may have different properties)
    struct Seed {
        uint64_t kmer;
        int pos;
        uint32_t key_idx, df_bytes;
        bool common;
    };

    // Candidate transcript structure (from discovery phase)
    struct Cand {
        uint32_t tid;       // Transcript ID
        int32_t off_q;      // Offset quantized (position >> bin_shift)
        uint16_t votes;     // Number of seed votes
        uint8_t flags;      // Junction flags (if applicable)
    };

    // Constructor
    explicit ONTMapper(const IndexVX &ix, FailureStats &st);

    // ==================== CORE MAPPING INTERFACE (MapperBase) ====================

    /**
     * Map a single ONT read to the transcriptome
     * Will implement in Step 7
     */
    std::vector<Alignment> map_read_single(
        const std::string &qname,
        const std::string &seq) override;

    /**
     * Map paired-end reads (not typical for ONT, but included for interface)
     * ONT is single-molecule sequencing, so this will log an error
     */
    void map_paired(
        const std::string &r1,
        const std::string &r2,
        const std::string &out_path,
        uint64_t max_pairs = 0,
        bool do_rescue = false) override;

    /**
     * Map single-end ONT reads from FASTQ file
     * Will implement in Step 9
     */
    void map_single_end(
        const std::string &reads,
        const std::string &out_path,
        uint64_t max_reads = 0,
        bool do_rescue = false) override;

    // ==================== EPOCH MANAGEMENT ====================

    /**
     * Reset per-read state for next read
     * Will implement in Step 4
     */
    void next_read_epoch() override;

    // ==================== ONT-SPECIFIC METHODS ====================

    /**
     * Compute adaptive stride based on read length
     * - Short reads (50-500bp): stride=1-2 (dense seeding)
     * - Medium reads (500-5kb): stride=5-10 (sparse seeding)
     * - Long reads (5kb+): stride=20-50 (very sparse seeding)
     * Will implement in Step 4
     */
    int compute_adaptive_stride(int read_len) const;

    /**
     * Generate seeds for ONT read using adaptive strategy
     * Will implement in Step 4
     */
    std::vector<Seed> make_seeds_ont(const std::string &seq) const;

    /**
     * Discovery result structure
     */
    struct DiscoverRes {
        std::vector<Cand> shortlist;
        uint16_t best_votes = 0;
        uint16_t best_runner = 0;
        bool had_any_seed = false;
        bool all_seeds_filtered = false;
    };


    // ==================== INHERITED INTERFACE IMPLEMENTATIONS ====================

    /**
     * Convert transcript coordinates to genomic coordinates
     * Can reuse Illumina's implementation (same reference genome)
     */
    GenomicPos transcript_to_genomic(uint32_t tid, int tx_pos, int read_len) const override;

    /**
     * Calculate MAPQ for ONT alignment
     * Will customize for long-read characteristics in Step 8
     */
    uint8_t calculate_mapq(const Alignment *best, const std::vector<Alignment> &all_hits, int read_len) const override;

    /**
     * Get minimum insert size (not relevant for ONT single-molecule)
     */
    int get_min_ins() const override { return min_ins; }

    /**
     * Get maximum insert size (not relevant for ONT single-molecule)
     */
    int get_max_ins() const override { return max_ins; }

private:
    // ==================== ALGORITHM COMPONENTS ====================

    // Hamming distance aligner (reused from Illumina, extended for ONT)
    Hamming2B ham;

    // ==================== ONT-SPECIFIC PARAMETERS ====================

    // Insert size parameters (not used for ONT, but required by interface)
    int min_ins = 100;
    int max_ins = 25000;  // Much larger than Illumina to account for long reads

    // Seeding parameters
    int min_votes = 2;      // Relaxed from Illumina's 3 (fewer seeds expected)
    int lead_margin = 1;    // Relaxed from Illumina's 2
    int base_stride = 10;   // Base stride for adaptive calculation
    int bin_shift = 3;      // Binning shift for position quantization

    // Error tolerance
    double max_error_rate = 0.20;  // 20% error tolerance for verification (vs Illumina's 10%)

    // Output options
    bool output_nm_md = false;  // Output NM and MD tags in SAM

    // ==================== DISCOVERY DATA STRUCTURES ====================

    /**
     * Heavy hitter structure for offset tracking
     * Tracks up to 4 most common offsets per transcript
     */
    struct HH4 {
        int32_t off_q[4];    // Quantized offsets
        uint16_t cnt[4];     // Vote counts
        uint8_t flg[4];      // Junction flags
        uint32_t gen;        // Generation counter

        HH4() : gen(0) {
            for (int i = 0; i < 4; i++) {
                off_q[i] = INT32_MIN;
                cnt[i] = 0;
                flg[i] = 0;
            }
        }
    };

    /**
     * Per-transcript best vote tracking
     */
    struct TBest {
        uint16_t best;          // Best vote count
        uint16_t second;        // Second-best vote count
        int32_t best_off_q;     // Best offset
        uint32_t gen;           // Generation counter

        TBest() : best(0), second(0), best_off_q(INT32_MIN), gen(0) {}
    };

    // Discovery state (mutable for const methods)
    mutable std::vector<HH4> hh;              // Heavy hitters per transcript
    mutable std::vector<uint32_t> hh_touched; // Touched transcripts
    mutable std::vector<TBest> per_tid;       // Best per transcript
    mutable uint32_t per_tid_gen;             // Generation counter

    // Sieve structures (for large candidate filtering)
    mutable std::vector<uint16_t> sieve_cnt;  // Seed counts per tid
    mutable std::vector<uint32_t> sieve_gen;  // Generation per tid
    mutable std::vector<uint32_t> tid_list;   // Candidate list

    // ==================== HELPER METHODS ====================

    /**
     * Generate seeds using simple adaptive stride
     * @param seq Read sequence
     * @param stride Stride between k-mers
     * @return Vector of seeds
     */
    std::vector<Seed> generate_seeds_simple(const std::string &seq, int stride) const;

    /**
     * Generate seeds using windowed strategy for ultra-long reads
     * @param seq Read sequence
     * @param base_stride Base stride value (will be modulated by windows)
     * @return Vector of seeds
     */
    std::vector<Seed> generate_seeds_windowed(const std::string &seq, int base_stride) const;

    /**
     * Discovery helper methods
     */
    uint16_t bump_hh(uint32_t tid, int32_t off_q, uint8_t flags);
    void on_vote(uint32_t tid, int32_t off_q, uint16_t cnt);

    /**
     * Discover candidate transcripts from ONT read
     * @param seq Read sequence (forward or reverse)
     * @param is_rev Is this the reverse strand?
     * @return Discovery result with shortlist of candidates
     */
    DiscoverRes discover_ont(const std::string &seq, bool is_rev);

    /**
     * Verify candidate alignment using extended Hamming distance
     * @param read Read sequence (forward or reverse)
     * @param read_2bit Read in 2-bit encoding
     * @param c Candidate from discovery
     * @param strand_is_rev Is this the reverse strand?
     * @return Alignment if valid, nullopt otherwise
     */
    std::optional<Alignment> verify_one_ont(
        const std::string &read,
        const std::vector<uint8_t> &read_2bit,
        const Cand &c,
        bool strand_is_rev) const;

    /**
     * Choose number of candidates to verify based on vote distribution
     * @param candidates List of candidates
     * @param best_votes Vote count of best candidate
     * @param best_runner Vote count of second-best candidate
     * @return Number of candidates to verify
     */
    int choose_verify_budget_ont(
        const std::vector<Cand> &candidates,
        uint16_t best_votes,
        uint16_t best_runner) const;
};

} // namespace rnamapper
