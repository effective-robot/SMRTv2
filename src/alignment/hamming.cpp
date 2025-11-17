// hamming.cpp - Hamming distance implementation
#include "hamming.h"
#include "../core/dna_utils.h"
#include <cstdint>

namespace rnamapper {

Hamming2B::Hamming2B() {
    for (int x = 0; x < 256; x++) {
        int m = 0;
        for (int g = 0; g < 4; g++) {
            int sh = 6 - 2 * g;
            if (((x >> sh) & 0x3) != 0) m++;
        }
        byte_pop4[x] = (uint8_t)m;
    }
}

int Hamming2B::run(const std::vector<uint8_t> &read_packed,
                   const uint8_t *target_packed,
                   uint32_t target_start,
                   int read_len,
                   int max_mm) const {
    int mm = 0;
    int full = read_len >> 2, rem = read_len & 3;
    uint32_t tbyte = target_start >> 2;
    int phase = target_start & 3;

    if (phase == 0) {
        // Aligned case - use SIMD-friendly 64-bit operations
        int full_bytes = full;
        int n64 = full_bytes >> 3;
        const uint64_t *rp = reinterpret_cast<const uint64_t *>(&read_packed[0]);
        const uint64_t *tp = reinterpret_cast<const uint64_t *>(&target_packed[tbyte]);
        for (int i = 0; i < n64; i++) {
            uint64_t x = rp[i] ^ tp[i];
            uint64_t y = (x | (x >> 1)) & 0x5555555555555555ULL;
            mm += __builtin_popcountll(y);
            if (mm > max_mm) return mm;
        }
        int consumed = n64 << 3;
        for (int i = consumed; i < full_bytes; i++) {
            uint8_t x = read_packed[i] ^ target_packed[tbyte + i];
            uint8_t y = (uint8_t)((x | (x >> 1)) & 0x55);
            mm += __builtin_popcount((unsigned)y);
            if (mm > max_mm) return mm;
        }
    } else {
        // Unaligned case
        int sh = phase * 2;
        for (int i = 0; i < full; i++) {
            uint8_t t = (uint8_t)((target_packed[tbyte + i] << sh) |
                                  (target_packed[tbyte + i + 1] >> (8 - sh)));
            uint8_t x = read_packed[i] ^ t;
            mm += byte_pop4[x];
            if (mm > max_mm) return mm;
        }
    }

    // Handle remainder bases
    for (int r = 0; r < rem; r++) {
        uint32_t ridx = (full << 2) + r;
        uint8_t rb = (read_packed[ridx >> 2] >> (6 - 2 * (ridx & 3))) & 0x3;
        uint8_t tb = get2b(target_packed, target_start + ridx);
        if (rb != tb) {
            mm++;
            if (mm > max_mm) return mm;
        }
    }
    return mm;
}

} // namespace rnamapper
