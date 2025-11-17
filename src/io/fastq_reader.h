// fastq_reader.h - FASTQ file reading with gzip support
#pragma once

#include <string>
#include <vector>
#include <cstdio>
#include <zlib.h>
#include "../core/types.h"

namespace rnamapper {

// Low-level line reader supporting both plain and gzipped files
class LineReader {
private:
    bool gz;
    FILE *fp;
    gzFile gzfp;
    std::vector<char> buf;

public:
    explicit LineReader(const std::string &path, size_t bufsz = 1 << 18);
    ~LineReader();

    // Delete copy constructor and assignment
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;

    bool getline(std::string &out);
};

// FASTQ format reader
class FastqReader {
private:
    LineReader lr;

public:
    explicit FastqReader(const std::string &path);

    bool next(std::string &id, std::string &seq, std::string &plus, std::string &qual);
    size_t next_batch(std::vector<ReadRecord> &batch, size_t max_reads);
};

} // namespace rnamapper
