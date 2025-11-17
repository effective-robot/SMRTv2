// index.h - Main index structure with k-mer lookup
#pragma once

#include "index_types.h"
#include <vector>
#include <string>
#include <utility>
#include <unordered_map>
#include <cstdint>
#include <functional>

namespace rnamapper {

class IndexVX {
public:
    // Index parameters
    uint32_t k = 31, flank = 50, prefix_bits = 20;
    uint32_t n_tx = 0, n_jx = 0, n_targets = 0;
    uint32_t n_buckets = 0, n_keys = 0;
    bool is_v12_plus = false;
    bool has_bloom = false;
    bool uses_canonical = false;

    // Metadata
    std::vector<TranscriptMeta> tx;
    std::vector<JunctionMeta> jx;
    std::vector<std::pair<std::string, uint32_t>> contigs;
    std::unordered_map<std::string, uint32_t> contig_len;
    std::unordered_map<std::string, uint32_t> tx_index_by_id;

    // Raw data pointers (non-owning views into blob)
    const uint32_t *bucket_cum = nullptr;
    const uint8_t *key_suffixes = nullptr;
    const uint32_t *post_offsets = nullptr;
    const uint8_t *common_bits = nullptr;
    const uint8_t *postings = nullptr;
    const uint8_t *seq_blob = nullptr;
    uint64_t seq_blob_size = 0;

    // Storage (owned data)
    std::vector<uint8_t> blob;
    std::vector<uint32_t> bucket_cum_storage;
    std::vector<uint8_t> suffix_storage;
    std::vector<uint32_t> post_off_storage;
    std::vector<uint8_t> common_storage;
    std::vector<uint64_t> bloom_bits;
    uint64_t bloom_n_bits = 0;
    int bloom_k_hashes = 3;

    // Load index from file
    void load(const std::string &path);

    // K-mer lookup functions (inline for performance)
    inline uint32_t prefix_of(uint64_t kmer) const {
        return (uint32_t)(kmer >> (62 - prefix_bits));
    }

    inline std::pair<uint32_t, uint32_t> bucket_range(uint32_t pre) const {
        return {bucket_cum[pre], bucket_cum[pre + 1]};
    }

    inline bool is_common(uint32_t key_idx) const {
        if (!common_bits) return false;
        return (common_bits[key_idx >> 3] >> (key_idx & 7)) & 1u;
    }

    inline uint64_t suffix_value_at(uint32_t idx) const {
        const uint8_t *a = key_suffixes + (size_t)idx * 6u;
        return (uint64_t)a[0] | ((uint64_t)a[1] << 8) | ((uint64_t)a[2] << 16) |
               ((uint64_t)a[3] << 24) | ((uint64_t)a[4] << 32) | ((uint64_t)a[5] << 40);
    }

    inline uint64_t suffix_of_kmer(uint64_t kmer) const {
        int bits = 62 - (int)prefix_bits;
        if (bits <= 0) return 0ULL;
        return kmer & ((bits == 64) ? ~0ULL : ((1ULL << bits) - 1ULL));
    }

    inline uint32_t postings_span_bytes(uint32_t key_idx) const {
        return post_offsets[key_idx + 1] - post_offsets[key_idx];
    }

    // Bloom filter functions (inline for performance)
    inline uint64_t bloom_hash1(uint64_t key) const {
        key ^= key >> 33;
        key *= 0xff51afd7ed558ccdULL;
        key ^= key >> 33;
        key *= 0xc4ceb9fe1a85ec53ULL;
        key ^= key >> 33;
        return key;
    }

    inline uint64_t bloom_hash2(uint64_t key) const {
        key = ~key + (key << 21);
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8);
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >> 28);
        key = key + (key << 31);
        return key;
    }

    inline bool bloom_may_contain(uint64_t kmer) const {
        if (!has_bloom) return true;
        uint64_t h1 = bloom_hash1(kmer), h2 = bloom_hash2(kmer);
        for (int i = 0; i < bloom_k_hashes; ++i) {
            uint64_t h = (h1 + i * h2) % bloom_n_bits;
            if (!(bloom_bits[h >> 6] & (1ULL << (h & 63))))
                return false;
        }
        return true;
    }

    // Find key in index
    std::pair<bool, uint32_t> find_key(uint64_t kmer) const;

    // Iterate over postings for a key (inline template for performance)
    template <typename CB>
    inline void for_each_posting_breakable(uint32_t key_idx, CB cb, bool do_prefetch) const {
        const uint32_t off0 = post_offsets[key_idx];
        const uint32_t off1 = post_offsets[key_idx + 1];
        const uint8_t *p = postings + off0;
        const uint8_t *e = postings + off1;

        auto zigzag_decode = [](uint64_t n) -> int64_t {
            return (int64_t)((n >> 1) ^ (-(int64_t)(n & 1)));
        };
        auto rd = [&](uint64_t &x) -> bool {
            x = 0;
            uint32_t s = 0;
            while (p < e) {
                if (do_prefetch && ((uintptr_t)p & 63u) == 0)
                    __builtin_prefetch(p + 128, 0, 1);
                uint8_t b = *p++;
                x |= (uint64_t)(b & 0x7FUL) << s;
                if (!(b & 0x80u)) return true;
                s += 7;
            }
            return false;
        };

        if (p < e && *p == 0xFF) {
            p++;
            if (e - p < 6) return;
            uint32_t tid = ((uint32_t)p[0]) | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16);
            uint32_t pos = ((uint32_t)p[3]) | ((uint32_t)p[4] << 8) | ((uint32_t)p[5] << 16);
            cb(tid, pos);
            return;
        }

        if (p >= e) return;
        uint64_t G = 0;
        if (!rd(G)) return;
        uint64_t prev_tid = 0;
        for (uint64_t g = 0; g < G && p < e; ++g) {
            uint64_t dt = 0, M = 0;
            if (!rd(dt)) return;
            uint32_t tid = (uint32_t)(prev_tid + dt);
            prev_tid = tid;
            if (!rd(M)) return;
            if (M == 0) continue;
            uint64_t first_pos_val = 0;
            if (!rd(first_pos_val)) return;
            uint32_t pos = (uint32_t)first_pos_val;
            if (!cb(tid, pos)) return;
            if (M > 1) {
                int64_t prev_delta = 0;
                for (uint64_t i = 1; i < M && p < e; ++i) {
                    uint64_t dod_enc = 0;
                    if (!rd(dod_enc)) return;
                    int64_t dod = zigzag_decode(dod_enc);
                    int64_t delta = prev_delta + dod;
                    pos = (uint32_t)((int64_t)pos + delta);
                    if (!cb(tid, pos)) return;
                    prev_delta = delta;
                }
            }
        }
    }

    // Helper functions
    inline bool is_jx_tid(uint32_t tid) const { return tid >= n_tx; }
    inline uint32_t jx_index(uint32_t tid) const { return tid - n_tx; }

private:
    // Binary reading helpers
    static inline uint32_t rd32(const uint8_t *&p, const uint8_t *end);
    static inline uint64_t rd64(const uint8_t *&p, const uint8_t *end);
    static inline std::string rd_str(const uint8_t *&p, const uint8_t *end);
};

} // namespace rnamapper
