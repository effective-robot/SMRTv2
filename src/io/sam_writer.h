// sam_writer.h - SAM format output
#pragma once

#include <string>
#include "../core/types.h"
#include "output_buffer.h"

// Forward declarations
namespace rnamapper {
class IndexVX;
class Mapper;
}

namespace rnamapper {

class SAMWriter {
private:
    const IndexVX &index;
    const Mapper &mapper;
    bool output_nm_md;

    // Helper functions
    static std::string normalize_qname(const std::string &s);
    void append_nm_md_tags(std::string &line, const Alignment *aln, const ReadRecord &rec) const;

public:
    SAMWriter(const IndexVX &ix, const Mapper &mp, bool nm_md = false);

    // Write SAM header with contigs and metadata
    void write_header(OutputBuffer &outbuf) const;

    // Write a single SAM record
    void write_record(OutputBuffer &buf, const ReadRecord &rec, const Alignment *aln,
                     bool is_paired, bool is_r2, const Alignment *mate_aln,
                     const ReadRecord *mate_rec, const std::vector<Alignment> *all_hits = nullptr);
};

} // namespace rnamapper
