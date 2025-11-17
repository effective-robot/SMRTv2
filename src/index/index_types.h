// index_types.h - Index metadata structures
#pragma once

#include <string>
#include <vector>
#include <utility>
#include <cstdint>

namespace rnamapper {

// Transcript metadata
struct TranscriptMeta {
    std::string id, chrom;
    char strand = '+';
    uint32_t len_bp = 0;
    uint32_t bytes = 0;
    const uint8_t *seq2b = nullptr;
    std::vector<std::pair<uint32_t, uint32_t>> genomic_exons;
};

// Junction metadata
struct JunctionMeta {
    std::string id, gene_id, chrom;
    char strand = '+';
    uint32_t donor_pos = 0, acceptor_pos = 0;
    uint32_t len_bp = 0;
    uint32_t bytes = 0;
    const uint8_t *seq2b = nullptr;
};

} // namespace rnamapper
