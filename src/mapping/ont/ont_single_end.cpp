// ont_single_end.cpp - ONT FASTQ file processing
#include "ont_mapper.h"
#include "../../io/fastq_reader.h"
#include "../../io/output_buffer.h"
#include "../../io/sam_writer.h"
#include <cstdio>

namespace rnamapper {

void ONTMapper::map_single_end(
    const std::string &reads_path,
    const std::string &out_path,
    uint64_t max_reads,
    bool do_rescue)
{
    (void)do_rescue;  // Not implemented for ONT

    // Open FASTQ reader
    FastqReader fq(reads_path);

    // Open SAM writer
    OutputBuffer outbuf(out_path);
    SAMWriter sam_writer(IX, *this, output_nm_md);
    sam_writer.write_header(outbuf);

    fprintf(stderr, "[ONTMapper] Processing FASTQ: %s\n", reads_path.c_str());
    fprintf(stderr, "[ONTMapper] Output: %s\n", out_path.c_str());
    fprintf(stderr, "[ONTMapper] Max reads: %s\n", max_reads ? std::to_string(max_reads).c_str() : "unlimited");

    // Batch processing
    const size_t BATCH_SIZE = 10000;  // Smaller batches for long reads
    std::vector<ReadRecord> batch;
    uint64_t total_processed = 0;

    while (true) {
        size_t n = fq.next_batch(batch, BATCH_SIZE);
        if (n == 0) break;

        for (size_t i = 0; i < n; ++i) {
            total_processed++;

            // Check max_reads limit
            if (max_reads && total_processed > max_reads) {
                goto finish;
            }

            const auto &rec = batch[i];

            // Map the read
            auto hits = map_read_single(rec.id, rec.seq);

            // Get best alignment (or nullptr if unmapped)
            Alignment *best = hits.empty() ? nullptr : &hits[0];

            // Write SAM record
            sam_writer.write_record(outbuf, rec, best,
                                   false,      // is_paired
                                   false,      // is_r2
                                   nullptr,    // mate_aln
                                   nullptr,    // mate_rec
                                   &hits);     // all_hits for MAPQ

            // Progress reporting
            if (total_processed % 10000 == 0) {
                fprintf(stderr, "[ONTMapper] Processed: %lu reads\n", total_processed);
            }
        }
    }

finish:
    fprintf(stderr, "\n");
    fprintf(stderr, "======================================================\n");
    fprintf(stderr, "ONT MAPPING SUMMARY\n");
    fprintf(stderr, "======================================================\n");
    fprintf(stderr, "Total reads processed:    %lu\n", total_processed);
    fprintf(stderr, "Reads mapped:             %lu (%.2f%%)\n",
            ST.ok_mapped_reads,
            100.0 * ST.ok_mapped_reads / std::max(1UL, total_processed));
    fprintf(stderr, "Reads with shortlist:     %lu\n", ST.reads_with_shortlist);
    fprintf(stderr, "Reads with verified hit:  %lu\n", ST.reads_with_verified_hit);
    fprintf(stderr, "Seeds filtered (Bloom):   %lu\n", ST.f_all_seeds_filtered);
    fprintf(stderr, "Discovery failed:         %lu\n", ST.f_discover_no_shortlist);
    fprintf(stderr, "Verification failed:      %lu\n", ST.f_verify_failed);
    fprintf(stderr, "======================================================\n");
}

} // namespace rnamapper
