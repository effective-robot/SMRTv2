// mapper_genomic.cpp - Transcript to genomic coordinate projection
#include "mapper.h"
#include <algorithm>

namespace rnamapper {

std::string Mapper::strip_isoform(const std::string &s) {
    auto d = s.find('.');
    return d == std::string::npos ? s : s.substr(0, d);
}

bool Mapper::tid_compatible_same_gene(uint32_t a, uint32_t b) const {
    if (a == b) return true;
    bool aj = IX.is_jx_tid(a), bj = IX.is_jx_tid(b);
    if (!aj && !bj) {
        auto g1 = strip_isoform(IX.tx[a].id);
        auto g2 = strip_isoform(IX.tx[b].id);
        return g1 == g2;
    }
    if (aj && !bj) {
        auto gT = strip_isoform(IX.tx[b].id);
        auto gJ = IX.jx[IX.jx_index(a)].gene_id;
        return gT == gJ;
    }
    if (!aj && bj) {
        auto gT = strip_isoform(IX.tx[a].id);
        auto gJ = IX.jx[IX.jx_index(b)].gene_id;
        return gT == gJ;
    }
    const auto &J1 = IX.jx[IX.jx_index(a)];
    const auto &J2 = IX.jx[IX.jx_index(b)];
    return J1.gene_id == J2.gene_id;
}

GenomicPos Mapper::transcript_to_genomic(uint32_t tid, int tx_pos, int read_len) const {
    GenomicPos gp{"", 0, '+', "", 0, false};
    if (read_len <= 0) return gp;

    auto build_blocks = [&](const std::vector<std::pair<uint32_t, uint32_t>> &tx_exons,
                           char strand,
                           int tx_pos0,
                           int len) -> std::vector<std::pair<uint32_t, uint32_t>> {
        std::vector<std::pair<uint32_t, uint32_t>> blocks;
        int tx_cursor = 0;
        int remaining = len;
        for (size_t i = 0; i < tx_exons.size() && remaining > 0; ++i) {
            uint32_t ex_s = tx_exons[i].first;
            uint32_t ex_e = tx_exons[i].second;
            int ex_len = (int)(ex_e - ex_s);
            if (tx_pos0 >= tx_cursor + ex_len) {
                tx_cursor += ex_len;
                continue;
            }
            int off_in_ex = std::max(0, tx_pos0 - tx_cursor);
            int avail = ex_len - off_in_ex;
            int take = std::min(remaining, avail);
            uint32_t block_start0;
            if (strand == '+') {
                block_start0 = ex_s + (uint32_t)off_in_ex;
            } else {
                block_start0 = (ex_e - (uint32_t)off_in_ex) - (uint32_t)take;
            }
            blocks.push_back({block_start0, (uint32_t)take});
            remaining -= take;
            tx_pos0 += take;
            tx_cursor += ex_len;
        }
        std::sort(blocks.begin(), blocks.end(),
                 [](const auto &a, const auto &b) { return a.first < b.first; });
        return blocks;
    };

    std::vector<std::pair<uint32_t, uint32_t>> tx_exons;
    char ref_strand = '+';
    std::string rname;

    if (IX.is_jx_tid(tid)) {
        const auto &J = IX.jx[IX.jx_index(tid)];
        rname = J.chrom;
        ref_strand = J.strand;
        uint32_t left_s = (uint32_t)(J.donor_pos) - IX.flank;
        uint32_t left_e = (uint32_t)(J.donor_pos);
        uint32_t right_s = (uint32_t)(J.acceptor_pos - 1);
        uint32_t right_e = right_s + IX.flank;
        if (ref_strand == '+') {
            tx_exons = {{left_s, left_e}, {right_s, right_e}};
        } else {
            tx_exons = {{right_s, right_e}, {left_s, left_e}};
        }
    } else {
        const auto &T = IX.tx[tid];
        rname = T.chrom;
        ref_strand = T.strand;
        if (T.genomic_exons.empty()) return gp;
        if (ref_strand == '+') {
            tx_exons.assign(T.genomic_exons.begin(), T.genomic_exons.end());
        } else {
            tx_exons.reserve(T.genomic_exons.size());
            for (auto it = T.genomic_exons.rbegin(); it != T.genomic_exons.rend(); ++it)
                tx_exons.push_back(*it);
        }
    }

    auto blocks = build_blocks(tx_exons, ref_strand, tx_pos, read_len);
    if (blocks.empty()) return gp;

    gp.chrom = rname;
    gp.strand = ref_strand;
    gp.pos = blocks.front().first + 1;
    gp.cigar.clear();
    gp.ref_span = 0;

    for (size_t i = 0; i < blocks.size(); ++i) {
        if (i > 0) {
            uint32_t prev_end = blocks[i - 1].first + blocks[i - 1].second;
            uint32_t intron = blocks[i].first - prev_end;
            gp.cigar += std::to_string(intron) + "N";
            gp.ref_span += intron;
        }
        gp.cigar += std::to_string(blocks[i].second) + "M";
        gp.ref_span += blocks[i].second;
    }

    gp.valid = true;
    return gp;
}

} // namespace rnamapper
