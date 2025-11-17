// sam_writer.cpp - SAM format output implementation
#include "sam_writer.h"
#include "../index/index.h"
#include "../mapping/mapper_base.h"
#include "../core/dna_utils.h"
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace rnamapper {

SAMWriter::SAMWriter(const IndexVX &ix, const MapperBase &mp, bool nm_md)
    : index(ix), mapper(mp), output_nm_md(nm_md) {}

std::string SAMWriter::normalize_qname(const std::string &s) {
    size_t len = s.size();
    if (len >= 2 && s[len - 2] == '/' && (s[len - 1] == '1' || s[len - 1] == '2')) {
        return s.substr(0, len - 2);
    }
    return s;
}

void SAMWriter::append_nm_md_tags(std::string &line, const Alignment *aln, const ReadRecord &rec) const {
    if (!aln || !output_nm_md) return;

    const uint8_t *ref_seq = index.is_jx_tid(aln->tid)
                                 ? index.jx[index.jx_index(aln->tid)].seq2b
                                 : index.tx[aln->tid].seq2b;
    int pos0 = aln->pos;
    int L = (int)rec.seq.size();
    int nm = 0;
    int run = 0;
    std::string md;
    md.reserve(L * 2);
    const char *nt = "ACGT";
    for (int i = 0; i < L; ++i) {
        char read_base = rec.seq[i];
        if (read_base == 'a')
            read_base = 'A';
        else if (read_base == 'c')
            read_base = 'C';
        else if (read_base == 'g')
            read_base = 'G';
        else if (read_base == 't' || read_base == 'u' || read_base == 'U')
            read_base = 'T';
        uint8_t ref_2b = get2b(ref_seq, pos0 + i);
        char ref_base = nt[ref_2b];
        if (read_base == ref_base) {
            ++run;
        } else {
            md += std::to_string(run);
            md += ref_base;
            run = 0;
            ++nm;
        }
    }
    md += std::to_string(run);
    line += "\tNM:i:";
    line += std::to_string(nm);
    line += "\tMD:Z:";
    line += md;
}

void SAMWriter::write_header(OutputBuffer &outbuf) const {
    std::string header = "@HD\tVN:1.6\tSO:unsorted\n";
    if (!index.contigs.empty()) {
        for (const auto &[chr_name, chr_len] : index.contigs) {
            header += "@SQ\tSN:" + chr_name + "\tLN:" + std::to_string(chr_len) + "\n";
        }
    } else {
        std::unordered_map<std::string, uint32_t> chr_max_pos;
        for (const auto &t : index.tx) {
            if (!t.genomic_exons.empty()) {
                uint32_t max_end = 0;
                for (const auto &[s, e] : t.genomic_exons) {
                    if (e > max_end) max_end = e;
                }
                uint32_t &current_max = chr_max_pos[t.chrom];
                if (max_end > current_max) current_max = max_end;
            }
        }
        for (const auto &j : index.jx) {
            uint32_t jx_end = std::max(j.donor_pos, j.acceptor_pos) + index.flank + 100;
            uint32_t &current_max = chr_max_pos[j.chrom];
            if (jx_end > current_max) current_max = jx_end;
        }
        std::vector<std::pair<std::string, uint32_t>> sorted_chrs(chr_max_pos.begin(), chr_max_pos.end());
        std::sort(sorted_chrs.begin(), sorted_chrs.end());
        for (const auto &[chr_name, chr_len] : sorted_chrs) {
            header += "@SQ\tSN:" + chr_name + "\tLN:" + std::to_string(chr_len) + "\n";
        }
    }
    header += "@RG\tID:rg1\tSM:sample\tPL:ILLUMINA\n";
    header += "@PG\tID:mapper_v17_profile_detailed\tPN:mapper_v17_profile_detailed\tVN:17\n";
    outbuf.append(header);
}

void SAMWriter::write_record(OutputBuffer &buf, const ReadRecord &rec, const Alignment *aln,
                            bool is_paired, bool is_r2, const Alignment *mate_aln,
                            const ReadRecord *mate_rec, const std::vector<Alignment> *all_hits) {
    std::string line;
    line.reserve(output_nm_md ? 384 : 256);
    line += normalize_qname(rec.id) + "\t";

    bool is_mapped = (aln != nullptr);
    bool mate_mapped = (mate_aln != nullptr);

    bool is_reverse = false, mate_reverse = false;
    if (aln) is_reverse = aln->is_reverse;
    if (mate_aln) mate_reverse = mate_aln->is_reverse;

    std::string rname = "*", mrname = "*";
    uint32_t pos = 0, mpos = 0;
    std::string cigar = "*";
    int read_len = (int)rec.seq.size();
    GenomicPos gpos, mate_gpos;

    if (is_mapped) {
        gpos = mapper.transcript_to_genomic(aln->tid, aln->pos, read_len);
        if (gpos.valid) {
            rname = gpos.chrom;
            pos = gpos.pos;
            cigar = gpos.cigar;
        } else {
            is_mapped = false;
            const_cast<FailureStats&>(mapper.ST).f_genomic_conv++;
            const_cast<FailureStats&>(mapper.ST).reads_with_tx_hit_but_bad_genomic++;
        }
    }

    if (mate_mapped) {
        int mate_len = mate_rec ? (int)mate_rec->seq.size() : read_len;
        mate_gpos = mapper.transcript_to_genomic(mate_aln->tid, mate_aln->pos, mate_len);
        if (mate_gpos.valid) {
            mrname = mate_gpos.chrom;
            mpos = mate_gpos.pos;
        } else {
            mate_mapped = false;
            const_cast<FailureStats&>(mapper.ST).f_genomic_conv++;
            const_cast<FailureStats&>(mapper.ST).reads_with_tx_hit_but_bad_genomic++;
        }
    }

    bool proper_pair = false;
    if (is_paired && is_mapped && mate_mapped && rname == mrname) {
        bool fr_orientation = (!is_reverse && mate_reverse && pos <= mpos) ||
                              (is_reverse && !mate_reverse && mpos <= pos);
        if (fr_orientation) {
            int ins = (int)pos + read_len - (int)mpos;
            if (ins < 0) ins = -ins;
            if (ins >= mapper.get_min_ins() && ins <= mapper.get_max_ins()) {
                proper_pair = true;
            }
        }
    }

    uint16_t flag = 0;
    if (is_paired) flag |= 0x1;
    if (proper_pair) flag |= 0x2;
    if (!is_mapped) flag |= 0x4;
    if (is_paired && !mate_mapped) flag |= 0x8;
    if (is_reverse) flag |= 0x10;
    if (is_paired && mate_reverse) flag |= 0x20;
    if (is_r2) flag |= 0x80;
    else if (is_paired) flag |= 0x40;

    line += std::to_string(flag) + "\t";

    if (is_mapped) {
        line += rname + "\t";
        line += std::to_string(pos) + "\t";
    } else {
        line += "*\t0\t";
    }

    // MAPQ calculation
    uint8_t mapq = 0;
    if (is_mapped) {
        if (all_hits && !all_hits->empty()) {
            mapq = mapper.calculate_mapq(aln, *all_hits, read_len);
        } else {
            if (aln->mm == 0) {
                mapq = 60;
            } else {
                double err_rate = (double)aln->mm / std::max(1, read_len);
                if (err_rate < 0.01) err_rate = 0.01;
                int q = (int)(-10.0 * log10(err_rate));
                mapq = (uint8_t)std::min(60, std::max(0, q));
            }
            if (aln->is_ambiguous) {
                mapq = std::min(mapq, (uint8_t)10);
            }
        }
    }

    line += std::to_string((int)mapq) + "\t";
    line += cigar + "\t";

    if (is_paired && mate_mapped) {
        if (rname == mrname) {
            line += "=";
        } else {
            line += mrname;
        }
        line.push_back('\t');
        line += std::to_string(mpos);
        line.push_back('\t');

        int tlen = 0;
        if (is_mapped && rname == mrname && mate_rec) {
            const int r1_len = read_len;
            const int r2_len = (int)mate_rec->seq.size();
            const int end1 = (int)pos + r1_len - 1;
            const int end2 = (int)mpos + r2_len - 1;
            const int left = std::min((int)pos, (int)mpos);
            const int right = std::max(end1, end2);
            tlen = right - left + 1;
            if ((int)pos != left) tlen = -tlen;
        }
        line += std::to_string(tlen);
        line.push_back('\t');
    } else {
        line += "*\t0\t0\t";
    }

    line += rec.seq;
    line.push_back('\t');
    line += rec.qual;
    line.push_back('\t');
    line += "RG:Z:rg1";

    if (is_mapped) {
        line.push_back('\t');
        line += "XS:A:";
        line.push_back(
            index.is_jx_tid(aln->tid)
                ? index.jx[index.jx_index(aln->tid)].strand
                : index.tx[aln->tid].strand);
        line.push_back('\t');
        line += "NH:i:1";
    }

    append_nm_md_tags(line, aln, rec);

    line.push_back('\n');
    buf.append(line);
}

} // namespace rnamapper
