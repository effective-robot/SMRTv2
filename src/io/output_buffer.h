// output_buffer.h - Buffered file output
#pragma once

#include <string>
#include <vector>
#include <cstdio>

namespace rnamapper {

class OutputBuffer {
private:
    std::vector<char> buf;
    size_t capacity;
    FILE *file = nullptr;

public:
    explicit OutputBuffer(const std::string &path, size_t cap = 32 * 1024 * 1024);
    ~OutputBuffer();

    // Delete copy constructor and assignment
    OutputBuffer(const OutputBuffer&) = delete;
    OutputBuffer& operator=(const OutputBuffer&) = delete;

    void append(const std::string &s);
    void flush();
};

} // namespace rnamapper
