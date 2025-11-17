// main.cpp - Entry point for RNA-seq mapper
#include "core/types.h"
#include "index/index.h"
#include "mapping/illumina/illumina_mapper.h"
#include <cstring>
#include <string>
#include <cstdio>

using namespace rnamapper;

// Helper functions for parsing CLI arguments
static inline bool ends_with(const std::string &s, const std::string &suf) {
    return s.size() >= suf.size() && std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}

static inline uint32_t parse_u32(const char *s, uint32_t dflt) {
    if (!s) return dflt;
    try {
        return (uint32_t)std::stoul(std::string(s));
    } catch (...) {
        return dflt;
    }
}

static inline int32_t parse_i32(const char *s, int32_t dflt) {
    if (!s) return dflt;
    try {
        return (int32_t)std::stol(std::string(s));
    } catch (...) {
        return dflt;
    }
}

static inline uint64_t parse_u64(const char *s, uint64_t dflt) {
    if (!s) return dflt;
    try {
        return (uint64_t)std::stoull(std::string(s));
    } catch (...) {
        return dflt;
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  Paired: %s <index.bin> <R1.fq[.gz]> <R2.fq[.gz]> <out.sam> [options]\n", argv[0]);
        fprintf(stderr, "  Single: %s <index.bin> <reads.fq[.gz]> <out.sam> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  --max-pairs N        limit paired mapping to N pairs\n");
        fprintf(stderr, "  --max N              alias for --max-pairs\n");
        fprintf(stderr, "  --stride N           seed stride\n");
        fprintf(stderr, "  --df1-bytes N        early DF cutoff\n");
        fprintf(stderr, "  --df2-bytes N        late DF cutoff\n");
        fprintf(stderr, "  --min-votes N        seed votes to pass\n");
        fprintf(stderr, "  --lead-margin N      lead over runner\n");
        fprintf(stderr, "  --early-seeds N      min seeds before early exit\n");
        fprintf(stderr, "  --verify-topk N      max verifies\n");
        fprintf(stderr, "  --verify-topk-lo N   verify when clear\n");
        fprintf(stderr, "  --output-nm-md 1     emit NM/MD\n");
        fprintf(stderr, "  --rescue             enable rescue/reprocess tiers for unmapped reads\n");
        return 1;
    }

    auto looks_like_fastq = [](const std::string &s) -> bool {
        return ends_with(s, ".fq") || ends_with(s, ".fastq") ||
               ends_with(s, ".fq.gz") || ends_with(s, ".fastq.gz");
    };

    std::string idx_path = argv[1];
    std::string arg2 = argv[2];
    std::string arg3 = argv[3];

    bool paired = false;
    std::string reads1 = arg2, reads2, out_path;
    int opt_start = 0;

    if (argc >= 5 && looks_like_fastq(arg3)) {
        paired = true;
        reads2 = arg3;
        out_path = argv[4];
        opt_start = 5;
    } else {
        paired = false;
        out_path = arg3;
        opt_start = 4;
    }

    uint64_t max_pairs = 0;
    uint64_t max_reads = 0;
    int min_ins = 100, max_ins = 800, stride = 1, min_votes = 3, lead_margin = 2, early_seeds = 6;
    int verify_topk = 12, verify_topk_lo = 5, sieve_seeds = 10, bin_shift = 3;
    uint32_t post_budget = 0, df1_bytes = 128, df2_bytes = 512;
    bool prefetch = true, dedup_seeds = false, skip_rev_on_clear = true;
    bool output_nm_md = false;
    bool rescue_unmapped = false;

    for (int i = opt_start; i < argc; i++) {
        if (strcmp(argv[i], "--max-pairs") == 0 && i + 1 < argc) {
            max_pairs = parse_u64(argv[++i], 0);
        } else if (strcmp(argv[i], "--max") == 0 && i + 1 < argc) {
            max_pairs = parse_u64(argv[++i], 0);
        } else if (strcmp(argv[i], "--max-reads") == 0 && i + 1 < argc) {
            max_reads = parse_u64(argv[++i], 0);
        } else if (strcmp(argv[i], "--min-ins") == 0 && i + 1 < argc) {
            min_ins = parse_i32(argv[++i], 100);
        } else if (strcmp(argv[i], "--max-ins") == 0 && i + 1 < argc) {
            max_ins = parse_i32(argv[++i], 800);
        } else if (strcmp(argv[i], "--stride") == 0 && i + 1 < argc) {
            stride = parse_i32(argv[++i], 1);
            if (stride < 1) stride = 1;
        } else if (strcmp(argv[i], "--df1-bytes") == 0 && i + 1 < argc) {
            df1_bytes = parse_u32(argv[++i], 128);
        } else if (strcmp(argv[i], "--df2-bytes") == 0 && i + 1 < argc) {
            df2_bytes = parse_u32(argv[++i], 512);
        } else if (strcmp(argv[i], "--post-budget") == 0 && i + 1 < argc) {
            post_budget = parse_u32(argv[++i], 0);
        } else if (strcmp(argv[i], "--min-votes") == 0 && i + 1 < argc) {
            min_votes = parse_i32(argv[++i], 3);
        } else if (strcmp(argv[i], "--lead-margin") == 0 && i + 1 < argc) {
            lead_margin = parse_i32(argv[++i], 2);
        } else if (strcmp(argv[i], "--early-seeds") == 0 && i + 1 < argc) {
            early_seeds = parse_i32(argv[++i], 6);
        } else if (strcmp(argv[i], "--verify-topk") == 0 && i + 1 < argc) {
            verify_topk = parse_i32(argv[++i], 12);
        } else if (strcmp(argv[i], "--verify-topk-lo") == 0 && i + 1 < argc) {
            verify_topk_lo = parse_i32(argv[++i], 5);
        } else if (strcmp(argv[i], "--prefetch") == 0 && i + 1 < argc) {
            prefetch = parse_i32(argv[++i], 1) != 0;
        } else if (strcmp(argv[i], "--dedup-seeds") == 0 && i + 1 < argc) {
            dedup_seeds = parse_i32(argv[++i], 0) != 0;
        } else if (strcmp(argv[i], "--sieve-seeds") == 0 && i + 1 < argc) {
            sieve_seeds = parse_i32(argv[++i], 10);
        } else if (strcmp(argv[i], "--skip-rev-on-clear") == 0 && i + 1 < argc) {
            skip_rev_on_clear = parse_i32(argv[++i], 1) != 0;
        } else if (strcmp(argv[i], "--bin-shift") == 0 && i + 1 < argc) {
            bin_shift = parse_i32(argv[++i], 3);
        } else if (strcmp(argv[i], "--output-nm-md") == 0 && i + 1 < argc) {
            output_nm_md = parse_i32(argv[++i], 0) != 0;
        } else if (strcmp(argv[i], "--rescue") == 0) {
            rescue_unmapped = true;
        }
    }

    IndexVX IX;
    IX.load(idx_path);

    FailureStats ST;
    IlluminaMapper mp(IX, ST);
    mp.min_ins = min_ins;
    mp.max_ins = max_ins;
    mp.stride = stride;
    mp.min_votes = min_votes;
    mp.lead_margin = lead_margin;
    mp.min_seeds_before_early = early_seeds;
    mp.verify_topk = verify_topk;
    mp.verify_topk_lo = verify_topk_lo;
    mp.df1_bytes = df1_bytes;
    mp.df2_bytes = df2_bytes;
    mp.post_budget = post_budget;
    mp.prefetch = prefetch;
    mp.dedup_seeds = dedup_seeds;
    mp.sieve_seeds = sieve_seeds;
    mp.skip_rev_on_clear = skip_rev_on_clear;
    mp.bin_shift = bin_shift;
    mp.output_nm_md = output_nm_md;

    if (!paired) {
        mp.map_single_end(reads1, out_path, max_reads, rescue_unmapped);
    } else {
        mp.map_paired(reads1, reads2, out_path, max_pairs, rescue_unmapped);
    }

    return 0;
}
