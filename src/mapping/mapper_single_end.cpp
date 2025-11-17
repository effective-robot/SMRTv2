// mapper_single_end.cpp - Single-end read mapping loop
#include "mapper.h"
#include "../io/sam_writer.h"
#include "../io/fastq_reader.h"
#include "../io/output_buffer.h"

namespace rnamapper {

void Mapper::map_single_end(const std::string &reads, const std::string &out_path,
                           uint64_t max_reads, bool do_rescue) {
    FastqReader fq(reads);
    OutputBuffer outbuf(out_path);
    SAMWriter sam_writer(IX, *this, output_nm_md);
    sam_writer.write_header(outbuf);

    uint64_t total = 0;
    const size_t BATCH_SIZE = 100000;
    std::vector<ReadRecord> batch;
    std::vector<Alignment> temp_alns;
    temp_alns.reserve(8);

    while (true) {
        size_t n = fq.next_batch(batch, BATCH_SIZE);
        if (n == 0) break;
        for (size_t idx = 0; idx < n; ++idx) {
            total++;
            if (max_reads && total > max_reads) break;

            const auto &rec = batch[idx];
            auto H = map_read_single(rec.id, rec.seq);
            Alignment *best = H.empty() ? nullptr : &H[0];

            if (do_rescue && !best) {
                auto rres = attempt_reprocess_single(rec.id, rec.seq);
                if (rres.recovered) {
                    temp_alns.clear();
                    temp_alns.push_back(rres.aln);
                    best = &temp_alns[0];
                }
            }

            if (best) ST.ok_mapped_reads++;
            sam_writer.write_record(outbuf, rec, best, false, false, nullptr, nullptr, &H);
        }
        if (max_reads && total >= max_reads) break;
    }
}

} // namespace rnamapper
