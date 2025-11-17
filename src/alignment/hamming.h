// hamming.h - Fast Hamming distance for 2-bit encoded sequences
#pragma once

#include <vector>
#include <cstdint>

namespace rnamapper {

class Hamming2B {
private:
    uint8_t byte_pop4[256];

public:
    Hamming2B();

    // Compute Hamming distance with early termination
    // Returns mismatch count (may exceed max_mm if early termination triggered)
    int run(const std::vector<uint8_t> &read_packed,
            const uint8_t *target_packed,
            uint32_t target_start,
            int read_len,
            int max_mm) const;
};

} // namespace rnamapper
