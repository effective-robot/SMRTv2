// fastq_reader.cpp
#include "fastq_reader.h"
#include <cstring>
#include <cstdlib>

namespace rnamapper {

// Helper function to check file extension
static inline bool ends_with(const std::string &s, const std::string &suf) {
    return s.size() >= suf.size() && std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}

LineReader::LineReader(const std::string &path, size_t bufsz)
    : gz(ends_with(path, ".gz")), fp(nullptr), gzfp(nullptr), buf(bufsz) {
    if (gz) {
        gzfp = gzopen(path.c_str(), "rb");
        if (!gzfp) {
            fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
            exit(1);
        }
        gzbuffer(gzfp, (unsigned int)bufsz);
    } else {
        fp = fopen(path.c_str(), "rb");
        if (!fp) {
            fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
            exit(1);
        }
        setvbuf(fp, nullptr, _IOFBF, bufsz);
    }
}

bool LineReader::getline(std::string &out) {
    out.clear();
    if (gz) {
        for (;;) {
            char *r = gzgets(gzfp, buf.data(), (int)buf.size());
            if (!r) return !out.empty();
            size_t len = strlen(r);
            if (len && r[len - 1] == '\n') {
                out.append(r, len - 1);
                return true;
            }
            out.append(r, len);
            if (len < buf.size() - 1) return true;
        }
    } else {
        int ch;
        while ((ch = fgetc(fp)) != EOF) {
            if (ch == '\n') return true;
            out.push_back((char)ch);
        }
        return !out.empty();
    }
}

LineReader::~LineReader() {
    if (gz && gzfp) gzclose(gzfp);
    if (!gz && fp) fclose(fp);
}

FastqReader::FastqReader(const std::string &path) : lr(path) {}

bool FastqReader::next(std::string &id, std::string &seq, std::string &plus, std::string &qual) {
    if (!lr.getline(id)) return false;
    if (!lr.getline(seq)) return false;
    if (!lr.getline(plus)) return false;
    if (!lr.getline(qual)) return false;

    // Strip @ and trim at first space
    if (!id.empty() && id[0] == '@') id = id.substr(1);
    auto sp = id.find_first_of(" \t");
    if (sp != std::string::npos) id = id.substr(0, sp);

    return true;
}

size_t FastqReader::next_batch(std::vector<ReadRecord> &batch, size_t max_reads) {
    batch.clear();
    batch.reserve(max_reads);
    for (size_t i = 0; i < max_reads; ++i) {
        ReadRecord r;
        if (!next(r.id, r.seq, r.plus, r.qual)) break;
        batch.push_back(std::move(r));
    }
    return batch.size();
}

} // namespace rnamapper
