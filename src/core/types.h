// types.h - Core data structures and type definitions
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace rnamapper {

// ======================= STATISTICS ==========================
struct FailureStats {
    uint64_t reprocess_attempts = 0;
    uint64_t reprocess_recovered = 0;
    uint64_t reprocess_failed = 0;

    uint64_t total_reads = 0;
    uint64_t total_pairs = 0;
    uint64_t ok_mapped_reads = 0;
    uint64_t ok_pairs = 0;
    uint64_t f_len_lt_k = 0;
    uint64_t f_no_seeds = 0;
    uint64_t f_all_seeds_filtered = 0;
    uint64_t f_discover_no_shortlist = 0;
    uint64_t f_verify_failed = 0;
    uint64_t reverse_attempted = 0;
    uint64_t reverse_skipped_clear = 0;
    uint64_t f_genomic_conv = 0;
    uint64_t f_pair_both_unmapped = 0;
    uint64_t f_pair_one_unmapped = 0;
    uint64_t f_pair_incompatible = 0;
    uint64_t reads_with_shortlist = 0;
    uint64_t reads_with_verified_hit = 0;
    uint64_t reads_with_tx_hit_but_bad_genomic = 0;
};

// ======================= READ RECORD ==========================
struct ReadRecord {
    std::string id;
    std::string seq;
    std::string plus;
    std::string qual;
};

// ======================= ALIGNMENT ==========================
struct Alignment {
    uint32_t tid;
    int pos;
    int mm;
    int score;
    bool is_reverse;
    bool is_ambiguous = false;
};

// ======================= GENOMIC POSITION ====================
struct GenomicPos {
    std::string chrom;
    uint32_t pos;
    char strand;
    std::string cigar;
    uint32_t ref_span;
    bool valid;
};

} // namespace rnamapper
