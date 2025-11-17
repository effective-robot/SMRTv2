// dna_utils.h - DNA sequence manipulation utilities
#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <cctype>

namespace rnamapper {

// Convert to uppercase and U->T
inline std::string upper(const std::string &s) {
    std::string out(s);
    for (char &c : out) {
        unsigned char u = toupper((unsigned char)c);
        c = (u == 'U') ? 'T' : u;
    }
    return out;
}

// Reverse complement
inline std::string revcomp(const std::string &s) {
    auto rc = [](char c) -> char {
        switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
        }
    };
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i)
        out[i] = rc(s[s.size() - 1 - i]);
    return out;
}

// Nucleotide to 2-bit encoding
inline int nt2b(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

// Pack sequence into 2-bit representation (4 bases per byte)
inline std::vector<uint8_t> pack2(const std::string &seq) {
    std::vector<uint8_t> out((seq.size() + 3) / 4, 0);
    size_t o = 0;
    int sh = 6;
    for (char c : seq) {
        int b = nt2b(c);
        if (b < 0) b = 0;
        out[o] |= (uint8_t)(b << sh);
        sh -= 2;
        if (sh < 0) {
            sh = 6;
            ++o;
        }
    }
    return out;
}

// Extract 2-bit base at given index
inline uint8_t get2b(const uint8_t *p, uint32_t idx) {
    uint32_t byte = idx >> 2, sh = 6 - 2 * (idx & 3);
    return (p[byte] >> sh) & 0x3;
}

} // namespace rnamapper
