// output_buffer.cpp
#include "output_buffer.h"
#include <cstdio>
#include <cstdlib>

namespace rnamapper {

OutputBuffer::OutputBuffer(const std::string &path, size_t cap) {
    capacity = cap;
    buf.reserve(capacity);
    file = fopen(path.c_str(), "wb");
    if (!file) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        exit(1);
    }
    setvbuf(file, nullptr, _IOFBF, 1 << 20);
}

void OutputBuffer::append(const std::string &s) {
    if (buf.size() + s.size() > capacity)
        flush();
    buf.insert(buf.end(), s.begin(), s.end());
}

void OutputBuffer::flush() {
    if (!buf.empty() && file) {
        fwrite(buf.data(), 1, buf.size(), file);
        buf.clear();
    }
}

OutputBuffer::~OutputBuffer() {
    flush();
    if (file) {
        fclose(file);
        file = nullptr;
    }
}

} // namespace rnamapper
