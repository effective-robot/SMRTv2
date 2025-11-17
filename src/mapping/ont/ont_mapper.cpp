// ont_mapper.cpp - ONTMapper implementation (constructor and basic setup)
#include "ont_mapper.h"
#include "../../core/dna_utils.h"
#include <cstdio>

namespace rnamapper {

ONTMapper::ONTMapper(const IndexVX &ix, FailureStats &st)
    : MapperBase(ix, st)
{
    fprintf(stderr, "[ONTMapper] Initialized for ONT long-read mode\n");
    fprintf(stderr, "[ONTMapper] Max error rate: %.1f%%, Base stride: %d\n",
            max_error_rate * 100.0, base_stride);
    // TODO Step 4: Initialize ONT-specific data structures
    // - Candidate hash tables
    // - Seed buffers
    // - Alignment caches
}

// ==================== CORE MAPPING INTERFACE (STUBS) ====================

std::vector<Alignment> ONTMapper::map_read_single(const std::string &qname, const std::string &seq) {
    (void)qname;  // Unused for now
    (void)seq;    // Unused for now
    fprintf(stderr, "ERROR: ONTMapper::map_read_single() not yet implemented (Step 7)\n");
    return {};
}

void ONTMapper::map_paired(const std::string &r1, const std::string &r2,
                           const std::string &out_path, uint64_t max_pairs, bool do_rescue) {
    (void)r1; (void)r2; (void)out_path; (void)max_pairs; (void)do_rescue;
    fprintf(stderr, "ERROR: ONTMapper::map_paired() not supported for ONT (single-molecule sequencing)\n");
    fprintf(stderr, "       Use map_single_end() instead or remove --ont flag for Illumina paired-end mode\n");
}

void ONTMapper::map_single_end(const std::string &reads, const std::string &out_path,
                               uint64_t max_reads, bool do_rescue) {
    (void)reads; (void)out_path; (void)max_reads; (void)do_rescue;
    fprintf(stderr, "ERROR: ONTMapper::map_single_end() not yet implemented (Step 9)\n");
}

// ==================== EPOCH MANAGEMENT ====================

void ONTMapper::next_read_epoch() {
    // TODO Step 4: Reset per-read state
    // - Clear candidate hash tables
    // - Reset seed counters
    // - Increment epoch generation
}

// ==================== ONT-SPECIFIC METHODS (STUBS) ====================

int ONTMapper::compute_adaptive_stride(int read_len) const {
    (void)read_len;  // Unused for now
    // TODO Step 4: Implement adaptive stride calculation
    // Strategy:
    //   - read_len < 500: stride = 1-2 (dense)
    //   - read_len < 5000: stride = 5-10 (sparse)
    //   - read_len >= 5000: stride = 20-50 (very sparse)
    return base_stride;  // Placeholder
}

std::vector<ONTMapper::Seed> ONTMapper::make_seeds_ont(const std::string &seq) const {
    (void)seq;  // Unused for now
    // TODO Step 4: Implement ONT seeding strategy
    // Options:
    //   1. Sparse k-mer sampling with adaptive stride
    //   2. Minimizer-based seeding
    //   3. Hybrid approach
    return {};
}

ONTMapper::DiscoverRes ONTMapper::discover_ont(const std::string &seq) const {
    (void)seq;  // Unused for now
    // TODO Step 5: Implement ONT candidate discovery
    // - Use sparse seeds
    // - Relaxed voting thresholds
    // - Handle high error rate
    DiscoverRes res;
    res.had_any_seed = false;
    return res;
}

bool ONTMapper::verify_ont(const std::string &read, Alignment &aln) const {
    (void)read;  // Unused for now
    (void)aln;   // Unused for now
    // TODO Step 6: Implement ONT-specific verification
    // - Use alignment algorithm tolerant of high error rate
    // - Compute alignment score
    // - Set mismatch count
    return false;
}

// ==================== INHERITED INTERFACE IMPLEMENTATIONS ====================

GenomicPos ONTMapper::transcript_to_genomic(uint32_t tid, int tx_pos, int read_len) const {
    (void)tid; (void)tx_pos; (void)read_len;
    // TODO Step 8: Implement or reuse Illumina's genomic projection
    // For now, return invalid position
    return {"", 0, '+', "", 0, false};
}

uint8_t ONTMapper::calculate_mapq(const Alignment *best, const std::vector<Alignment> &all_hits, int read_len) const {
    (void)best; (void)all_hits; (void)read_len;
    // TODO Step 8: Implement ONT-specific MAPQ calculation
    // Consider:
    // - Long read length
    // - Higher error rate
    // - Different alignment score distributions
    return 0;  // Placeholder
}

} // namespace rnamapper
