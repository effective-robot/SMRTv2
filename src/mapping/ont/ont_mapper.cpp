// ont_mapper.cpp - ONTMapper implementation (constructor and basic setup)
#include "ont_mapper.h"
#include "../../core/dna_utils.h"
#include <cstdio>
#include <algorithm>
#include <cmath>

namespace rnamapper {

ONTMapper::ONTMapper(const IndexVX &ix, FailureStats &st)
    : MapperBase(ix, st), per_tid_gen(1)
{
    fprintf(stderr, "[ONTMapper] Initialized for ONT long-read mode\n");
    fprintf(stderr, "[ONTMapper] Max error rate: %.1f%%, Base stride: %d\n",
            max_error_rate * 100.0, base_stride);

    // Initialize discovery structures
    per_tid.resize(IX.n_tx + IX.n_jx);
    hh.resize(IX.n_tx + IX.n_jx);
    hh_touched.reserve(512);
    sieve_cnt.resize(IX.n_targets, 0);
    sieve_gen.resize(IX.n_targets, 0);
    tid_list.reserve(256);
}

// ==================== CORE MAPPING INTERFACE (STUBS) ====================

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

// ==================== ONT-SPECIFIC METHODS ====================

std::vector<ONTMapper::Seed> ONTMapper::make_seeds_ont(const std::string &seq) const {
    const int L = (int)seq.size();
    const int K = (int)IX.k;  // k=31 from index

    std::vector<Seed> seeds;
    if (L < K) {
        return seeds;  // Read too short
    }

    // Step 1: Compute adaptive stride for this read (implemented in ont_seeding.cpp)
    int adaptive_stride = compute_adaptive_stride(L);

    // Step 2: Determine if windowing is needed (ultra-long reads)
    bool use_windowed = (L > 10000);  // Step 5 threshold

    // Step 3: Generate seeds
    if (!use_windowed) {
        // Normal adaptive stride seeding
        seeds = generate_seeds_simple(seq, adaptive_stride);
    } else {
        // Windowed seeding for ultra-long reads (Step 5)
        seeds = generate_seeds_windowed(seq, adaptive_stride);
    }

    // Step 4: Sort seeds by DF (document frequency) for efficiency
    std::sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
        if (a.common != b.common) return a.common < b.common;
        if (a.df_bytes != b.df_bytes) return a.df_bytes < b.df_bytes;
        return a.pos < b.pos;
    });

    // Debug logging
    if (use_windowed) {
        fprintf(stderr, "[ONT Seeding] Ultra-long read: %d bp, windowed mode, seeds: %zu\n",
                L, seeds.size());
    } else if (L > 1000) {
        fprintf(stderr, "[ONT Seeding] Read: %d bp, stride: %d, seeds: %zu\n",
                L, adaptive_stride, seeds.size());
    }

    return seeds;
}

std::vector<ONTMapper::Seed> ONTMapper::generate_seeds_simple(
    const std::string &seq, int stride) const
{
    const int L = (int)seq.size();
    const int K = (int)IX.k;
    const uint64_t KMASK = (K == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * K)) - 1ULL);

    std::vector<Seed> seeds;
    seeds.reserve(std::max(0, (L - K + stride) / stride));

    uint64_t roll = 0;
    int good = 0;

    auto push_seed = [&](int start, uint64_t kmer) {
        // Canonicalize if index uses canonical k-mers
        if (IX.uses_canonical) {
            uint64_t rc = 0, kx = kmer;
            for (int i = 0; i < K; i++) {
                uint64_t base = kx & 3ULL;
                rc = (rc << 2) | (3ULL - base);
                kx >>= 2;
            }
            if (rc < kmer) kmer = rc;
        }

        // Bloom filter check
        if (!IX.bloom_may_contain(kmer)) return;

        // Index lookup
        auto [ok, idx] = IX.find_key(kmer);
        if (!ok) return;

        // Get document frequency
        uint32_t byte_span = IX.postings_span_bytes(idx);
        uint32_t df_bytes = (byte_span <= 6u) ? 6u : byte_span;

        seeds.push_back(Seed{kmer, start, idx, df_bytes, IX.is_common(idx)});
    };

    // Rolling hash to extract k-mers
    for (int i = 0; i < L; i++) {
        int b = nt2b(seq[i]);
        if (b < 0) {
            roll = 0;
            good = 0;
            continue;
        }
        roll = ((roll << 2) | (uint64_t)b) & KMASK;
        if (++good >= K) {
            int start = i - K + 1;
            if ((start % stride) == 0) {
                push_seed(start, roll);
            }
        }
    }

    return seeds;
}

std::vector<ONTMapper::Seed> ONTMapper::generate_seeds_windowed(
    const std::string &seq, int base_stride) const
{
    /**
     * Windowed seeding strategy for ultra-long reads (>10kb)
     *
     * Strategy:
     * - Divide read into windows (e.g., 2kb windows every 5kb)
     * - Sample DENSELY within first 2kb of each window
     * - Sample SPARSELY between windows
     *
     * Goal: Get good coverage without excessive seed count
     */

    const int L = (int)seq.size();
    const int K = (int)IX.k;
    const uint64_t KMASK = (K == 32) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (2 * K)) - 1ULL);

    const int window_size = 2000;      // Dense sampling in 2kb windows
    const int window_stride = 5000;    // Windows every 5kb

    std::vector<Seed> seeds;
    seeds.reserve(L / 50);  // Estimate

    auto should_sample = [&](int pos) -> bool {
        int pos_in_window = pos % window_stride;

        if (pos_in_window < window_size) {
            // Dense sampling in window (use base_stride)
            return (pos % base_stride) == 0;
        } else {
            // Sparse sampling between windows (3x stride)
            return (pos % (base_stride * 3)) == 0;
        }
    };

    auto push_seed = [&](int start, uint64_t kmer) {
        // Canonicalize if index uses canonical k-mers
        if (IX.uses_canonical) {
            uint64_t rc = 0, kx = kmer;
            for (int i = 0; i < K; i++) {
                uint64_t base = kx & 3ULL;
                rc = (rc << 2) | (3ULL - base);
                kx >>= 2;
            }
            if (rc < kmer) kmer = rc;
        }

        // Bloom filter check
        if (!IX.bloom_may_contain(kmer)) return;

        // Index lookup
        auto [ok, idx] = IX.find_key(kmer);
        if (!ok) return;

        // Get document frequency
        uint32_t byte_span = IX.postings_span_bytes(idx);
        uint32_t df_bytes = (byte_span <= 6u) ? 6u : byte_span;

        seeds.push_back(Seed{kmer, start, idx, df_bytes, IX.is_common(idx)});
    };

    // Generate seeds with windowing logic
    uint64_t roll = 0;
    int good = 0;

    for (int i = 0; i < L; i++) {
        int b = nt2b(seq[i]);
        if (b < 0) {
            roll = 0;
            good = 0;
            continue;
        }
        roll = ((roll << 2) | (uint64_t)b) & KMASK;
        if (++good >= K) {
            int start = i - K + 1;
            if (should_sample(start)) {
                push_seed(start, roll);
            }
        }
    }

    return seeds;
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
