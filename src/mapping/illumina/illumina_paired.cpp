// illumina_paired.cpp - Paired-end scoring and mate rescue
#include "illumina_mapper.h"
#include "../../io/sam_writer.h"
#include "../../io/fastq_reader.h"
#include "../../io/output_buffer.h"
#include <algorithm>
#include <cstdlib>
#include <set>

namespace rnamapper {

IlluminaMapper::PairScore IlluminaMapper::score_pair(const Alignment &a1, const Alignment &a2,
                                     const ReadRecord &rec1, const ReadRecord &rec2) const {
    PairScore ps;
    ps.r1_score = a1.score;
    ps.r2_score = a2.score;
    ps.concordance_bonus = 0;
    ps.insert_quality = 0;

    GenomicPos gp1 = transcript_to_genomic(a1.tid, a1.pos, (int)rec1.seq.size());
    GenomicPos gp2 = transcript_to_genomic(a2.tid, a2.pos, (int)rec2.seq.size());

    if (gp1.valid && gp2.valid) {
        if (gp1.chrom == gp2.chrom) {
            int ins = abs((int)gp1.pos - (int)gp2.pos);

            if (ins >= min_ins && ins <= max_ins) {
                ps.concordance_bonus = 15;

                int ideal_ins = 400;
                int ins_dev = abs(ins - ideal_ins);
                if (ins_dev < 50) {
                    ps.insert_quality = 20;
                } else if (ins_dev < 150) {
                    ps.insert_quality = 10;
                } else {
                    ps.insert_quality = 5;
                }
            } else if (ins <= max_ins * 2) {
                ps.concordance_bonus = 20;
            }
        } else {
            ps.concordance_bonus = -15;
        }
    }

    ps.total = ps.r1_score + ps.r2_score + ps.concordance_bonus + ps.insert_quality;
    return ps;
}

IlluminaMapper::BestPair IlluminaMapper::select_best_pair(std::vector<Alignment> &H1, std::vector<Alignment> &H2,
                                          const ReadRecord &rec1, const ReadRecord &rec2) const {
    BestPair result;
    result.best1 = nullptr;
    result.best2 = nullptr;
    result.found = false;
    result.score.total = INT_MIN;

    if (H1.empty() || H2.empty()) return result;

    const int MAX_CHECK = 10;
    int n1 = std::min<int>(MAX_CHECK, (int)H1.size());
    int n2 = std::min<int>(MAX_CHECK, (int)H2.size());

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            PairScore ps = score_pair(H1[i], H2[j], rec1, rec2);

            if (ps.total > result.score.total) {
                result.score = ps;
                result.best1 = &H1[i];
                result.best2 = &H2[j];
                result.found = true;
            }
        }
    }

    return result;
}

bool IlluminaMapper::attempt_mate_rescue(std::vector<Alignment> &H1, Alignment *&best1,
                                 std::vector<Alignment> &H2, Alignment *&best2,
                                 const ReadRecord &rec1, const ReadRecord &rec2) {
    if (!best1 || !best2) return false;

    GenomicPos gp1 = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
    GenomicPos gp2 = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());

    if (!gp1.valid || !gp2.valid) return false;

    if (gp1.chrom == gp2.chrom) {
        int ins = abs((int)gp1.pos - (int)gp2.pos);
        if (ins >= min_ins && ins <= max_ins) return false;
    }

    struct RescueOption {
        Alignment *new_best1;
        Alignment *new_best2;
        int insert_size;
        int score_penalty;
        double quality;
    };

    std::vector<RescueOption> options;

    // Option 1: Keep R1, try alternatives for R2
    for (size_t i = 0; i < H2.size() && i < 10; i++) {
        GenomicPos alt_gp2 = transcript_to_genomic(H2[i].tid, H2[i].pos, (int)rec2.seq.size());
        if (alt_gp2.valid && alt_gp2.chrom == gp1.chrom) {
            int ins = abs((int)gp1.pos - (int)alt_gp2.pos);
            if (ins <= max_ins * 2) {
                int penalty = best2->score - H2[i].score;
                int ideal_ins = 400;
                int ins_dev = abs(ins - ideal_ins);
                double quality = 100.0 - penalty - (ins_dev / 10.0);
                options.push_back({best1, &H2[i], ins, penalty, quality});
            }
        }
    }

    // Option 2: Keep R2, try alternatives for R1
    for (size_t i = 0; i < H1.size() && i < 10; i++) {
        GenomicPos alt_gp1 = transcript_to_genomic(H1[i].tid, H1[i].pos, (int)rec1.seq.size());
        if (alt_gp1.valid && alt_gp1.chrom == gp2.chrom) {
            int ins = abs((int)alt_gp1.pos - (int)gp2.pos);
            if (ins <= max_ins * 2) {
                int penalty = best1->score - H1[i].score;
                int ideal_ins = 400;
                int ins_dev = abs(ins - ideal_ins);
                double quality = 100.0 - penalty - (ins_dev / 10.0);
                options.push_back({&H1[i], best2, ins, penalty, quality});
            }
        }
    }

    if (options.empty()) return false;

    std::sort(options.begin(), options.end(),
             [](const RescueOption &a, const RescueOption &b) {
                 return a.quality > b.quality;
             });

    const auto &best_option = options[0];

    if (best_option.insert_size <= max_ins &&
        best_option.score_penalty <= 10 &&
        best_option.quality > 0) {
        best1 = best_option.new_best1;
        best2 = best_option.new_best2;
        return true;
    }

    return false;
}

void IlluminaMapper::map_paired(const std::string &r1, const std::string &r2,
                       const std::string &out_path, uint64_t max_pairs, bool do_rescue) {
    // Structure for tracking ambiguous pairs that need dense remapping
    struct AmbiguousPair {
        std::string qname;
        std::string r1_seq;
        std::string r2_seq;
        std::string r1_qual;
        std::string r2_qual;
    };

    std::vector<AmbiguousPair> g_ambiguous_pairs;
    std::set<std::string> g_skip_pairs;

    FastqReader fq1(r1), fq2(r2);
    OutputBuffer outbuf(out_path);
    SAMWriter sam_writer(IX, *this, output_nm_md);
    sam_writer.write_header(outbuf);

    const size_t BATCH_SIZE = 50000;
    std::vector<ReadRecord> batch1, batch2;
    uint64_t total = 0;

    std::vector<Alignment> temp_alns1, temp_alns2;
    temp_alns1.reserve(8);
    temp_alns2.reserve(8);

    while (true) {
        size_t n1 = fq1.next_batch(batch1, BATCH_SIZE);
        size_t n2 = fq2.next_batch(batch2, BATCH_SIZE);
        if (n1 == 0 || n2 == 0) break;
        if (n1 != n2) {
            fprintf(stderr, "ERROR: FASTQ mismatch\n");
            exit(1);
        }

        for (size_t idx = 0; idx < n1; ++idx) {
            total++;
            ST.total_pairs++;

            if (max_pairs && total > max_pairs) goto finish;

            const auto &rec1 = batch1[idx];
            const auto &rec2 = batch2[idx];

            auto H1 = map_read_single(rec1.id, rec1.seq);
            auto H2 = map_read_single(rec2.id, rec2.seq);

            Alignment *best1 = nullptr;
            Alignment *best2 = nullptr;

            if (!H1.empty() && !H2.empty()) {
                best1 = &H1[0];
                best2 = &H2[0];

                bool has_alternatives = (H1.size() > 1 || H2.size() > 1);

                if (has_alternatives) {
                    GenomicPos gp1 = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                    GenomicPos gp2 = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());

                    bool discordant = (gp1.valid && gp2.valid && gp1.chrom != gp2.chrom);
                    bool r1_close = (H1.size() > 1 && H1[0].score - H1[1].score <= 10);
                    bool r2_close = (H2.size() > 1 && H2[0].score - H2[1].score <= 10);

                    if (discordant || r1_close || r2_close) {
                        auto best_pair = select_best_pair(H1, H2, rec1, rec2);

                        if (best_pair.found) {
                            int independent_score = H1[0].score + H2[0].score;
                            int pair_score = best_pair.score.total - best_pair.score.concordance_bonus;

                            if (pair_score >= independent_score - 10) {
                                best1 = best_pair.best1;
                                best2 = best_pair.best2;
                            }
                        }
                    }
                }
            } else {
                best1 = H1.empty() ? nullptr : &H1[0];
                best2 = H2.empty() ? nullptr : &H2[0];
            }

            // Mate rescue
            if (best1 && best2) {
                GenomicPos gp1_pre = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                GenomicPos gp2_pre = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());
                if (gp1_pre.valid && gp2_pre.valid && gp1_pre.chrom != gp2_pre.chrom) {
                    attempt_mate_rescue(H1, best1, H2, best2, rec1, rec2);
                }
            }

            // CRITICAL: Detect ambiguous pairs for dense remapping
            // This is the missing logic that caused the regression
            bool pair_is_ambiguous = false;
            if (best1 && best2) {
                bool r1_ambig = best1->is_ambiguous;
                bool r2_ambig = best2->is_ambiguous;

                GenomicPos gp1_check = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                GenomicPos gp2_check = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());

                bool discordant = false;
                if (gp1_check.valid && gp2_check.valid) {
                    if (gp1_check.chrom != gp2_check.chrom &&
                        gp1_check.chrom != "*" && gp2_check.chrom != "*") {
                        discordant = true;
                    }
                }

                if (discordant && (r1_ambig || r2_ambig)) {
                    pair_is_ambiguous = true;
                }
            }

            if (pair_is_ambiguous) {
                g_ambiguous_pairs.push_back({rec1.id, rec1.seq, rec2.seq, rec1.qual, rec2.qual});
                g_skip_pairs.insert(rec1.id);
            }

            // Rescue unmapped reads
            if (do_rescue && !best1) {
                auto rres = attempt_reprocess_single(rec1.id, rec1.seq);
                if (rres.recovered) {
                    temp_alns1.clear();
                    temp_alns1.push_back(rres.aln);
                    best1 = &temp_alns1[0];
                }
            }

            if (do_rescue && !best2) {
                auto rres = attempt_reprocess_single(rec2.id, rec2.seq);
                if (rres.recovered) {
                    temp_alns2.clear();
                    temp_alns2.push_back(rres.aln);
                    best2 = &temp_alns2[0];
                }
            }

            bool conc = false;
            if (best1 && best2) {
                if (tid_compatible_same_gene(best1->tid, best2->tid)) {
                    if (best1->tid == best2->tid) {
                        int ins = (best1->pos + (int)rec1.seq.size()) - best2->pos;
                        if (ins < 0) ins = -ins;
                        if (ins >= min_ins && ins <= max_ins) conc = true;
                    } else {
                        conc = true;
                    }
                }
            }

            if (!best1 && !best2) {
                ST.f_pair_both_unmapped++;
            } else if (!best1 || !best2) {
                ST.f_pair_one_unmapped++;
            } else {
                if (conc) {
                    ST.ok_pairs++;
                    ST.ok_mapped_reads += 2;
                } else {
                    ST.f_pair_incompatible++;
                    ST.ok_mapped_reads += 2;
                }
            }

            // Skip output for ambiguous pairs - they will be remapped and written later
            if (g_skip_pairs.find(rec1.id) == g_skip_pairs.end()) {
                sam_writer.write_record(outbuf, rec1, best1, true, false, best2, &rec2, &H1);
                sam_writer.write_record(outbuf, rec2, best2, true, true, best1, &rec1, &H2);
            }
        }
    }

finish:
    // Dense remapping of ambiguous pairs with stride=1
    if (!g_ambiguous_pairs.empty()) {
        fprintf(stderr, "Detected %zu ambiguous pairs, remapping with stride=1...\n",
                g_ambiguous_pairs.size());

        // Save current settings
        int old_stride = stride;
        int old_min_votes = min_votes;
        int old_lead_margin = lead_margin;
        uint32_t old_df1 = df1_bytes;
        uint32_t old_df2 = df2_bytes;
        bool old_skip_common = skip_common;
        int old_verify_topk = verify_topk;

        // Apply dense remapping settings
        stride = 1;
        min_votes = 2;
        lead_margin = 3;
        df1_bytes = 256;
        df2_bytes = 1024;
        skip_common = false;
        verify_topk = 30;

        size_t remapped_count = 0;
        for (const auto &amb : g_ambiguous_pairs) {
            remapped_count++;
            if (remapped_count % 10 == 0) {
                fprintf(stderr, "  Remapping ambiguous pairs: %zu/%zu\n",
                        remapped_count, g_ambiguous_pairs.size());
            }

            ReadRecord rec1{amb.qname, amb.r1_seq, "+", amb.r1_qual};
            ReadRecord rec2{amb.qname, amb.r2_seq, "+", amb.r2_qual};

            auto H1 = map_read_single(amb.qname, amb.r1_seq);
            auto H2 = map_read_single(amb.qname, amb.r2_seq);

            Alignment *best1 = nullptr;
            Alignment *best2 = nullptr;

            if (!H1.empty() && !H2.empty()) {
                best1 = &H1[0];
                best2 = &H2[0];

                bool has_alternatives = (H1.size() > 1 || H2.size() > 1);

                if (has_alternatives) {
                    GenomicPos gp1 = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                    GenomicPos gp2 = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());

                    bool discordant = (gp1.valid && gp2.valid && gp1.chrom != gp2.chrom);
                    bool r1_close = (H1.size() > 1 && H1[0].score - H1[1].score <= 10);
                    bool r2_close = (H2.size() > 1 && H2[0].score - H2[1].score <= 10);

                    if (discordant || r1_close || r2_close) {
                        auto best_pair = select_best_pair(H1, H2, rec1, rec2);

                        if (best_pair.found) {
                            int independent_score = H1[0].score + H2[0].score;
                            int pair_score = best_pair.score.total - best_pair.score.concordance_bonus;

                            if (pair_score >= independent_score - 10) {
                                best1 = best_pair.best1;
                                best2 = best_pair.best2;
                            }
                        }
                    }
                }
            } else {
                best1 = H1.empty() ? nullptr : &H1[0];
                best2 = H2.empty() ? nullptr : &H2[0];
            }

            // Mate rescue for remapped pairs
            if (best1 && best2) {
                GenomicPos gp1_pre = transcript_to_genomic(best1->tid, best1->pos, (int)rec1.seq.size());
                GenomicPos gp2_pre = transcript_to_genomic(best2->tid, best2->pos, (int)rec2.seq.size());
                if (gp1_pre.valid && gp2_pre.valid && gp1_pre.chrom != gp2_pre.chrom) {
                    attempt_mate_rescue(H1, best1, H2, best2, rec1, rec2);
                }
            }

            // Write remapped results
            sam_writer.write_record(outbuf, rec1, best1, true, false, best2, &rec2, &H1);
            sam_writer.write_record(outbuf, rec2, best2, true, true, best1, &rec1, &H2);
        }

        fprintf(stderr, "Finished remapping %zu ambiguous pairs.\n", g_ambiguous_pairs.size());

        // Restore original settings
        stride = old_stride;
        min_votes = old_min_votes;
        lead_margin = old_lead_margin;
        df1_bytes = old_df1;
        df2_bytes = old_df2;
        skip_common = old_skip_common;
        verify_topk = old_verify_topk;
    }
}

} // namespace rnamapper
