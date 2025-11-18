// ont_seeding.cpp - ONT seed generation and adaptive stride selection
#include "ont_mapper.h"
#include <algorithm>

namespace rnamapper {

int ONTMapper::compute_adaptive_stride(int read_len) const {
    /**
     * Adaptive stride strategy for ONT reads
     *
     * Goals:
     * - Short reads (50-200bp): Dense seeding (~50-100 seeds)
     * - Medium reads (200bp-5kb): Moderate seeding (~80-120 seeds)
     * - Long reads (5kb-22kb): Sparse seeding (~60-100 seeds)
     *
     * Formula: stride = max(1, read_len / target_seeds)
     *
     * Rationale:
     * - ONT has 10-20% error rate, but longer reads provide more signal
     * - Maintain reasonable seed count regardless of read length
     * - Balance: enough seeds for discovery, not too many for speed
     */

    if (read_len < 200) {
        // Short reads: dense sampling (stride=1-2)
        // ~50-100 seeds for 100bp read
        return 1;
    }
    else if (read_len < 1000) {
        // Medium-short: target ~100 seeds
        // e.g., 500bp read → stride 5 → 100 seeds
        return std::max(1, read_len / 100);
    }
    else if (read_len < 5000) {
        // Medium-long: target ~80 seeds
        // e.g., 3000bp read → stride 37 → 81 seeds
        return std::max(1, read_len / 80);
    }
    else if (read_len < 10000) {
        // Long: target ~70 seeds
        // e.g., 7000bp read → stride 100 → 70 seeds
        return std::max(1, read_len / 70);
    }
    else {
        // Ultra-long: target ~60 seeds (with windowing, see Step 5)
        // e.g., 20000bp read → stride 333 → 60 seeds
        // Windowed mode will add density at specific regions
        return std::max(1, read_len / 60);
    }
}

} // namespace rnamapper
