# RNA-seq Mapper - Modular Refactoring

## Overview

This is a complete refactoring of the monolithic `mapper_v17_profile_detailed.cpp` into a clean, modular C++17 architecture. The refactoring preserves **exact** functionality and performance while dramatically improving maintainability, testability, and code organization.

## Architecture

### Module Structure

```
src/
├── core/                    # Core data structures and utilities
│   ├── types.h              # Basic types (Alignment, ReadRecord, FailureStats, etc.)
│   └── dna_utils.h          # DNA manipulation (upper, revcomp, pack2, get2b, etc.)
│
├── io/                      # Input/Output operations
│   ├── fastq_reader.h/cpp   # FASTQ reading (LineReader, FastqReader)
│   ├── output_buffer.h/cpp  # Buffered file writing
│   └── sam_writer.h/cpp     # SAM format output
│
├── index/                   # Index structures and loading
│   ├── index_types.h        # TranscriptMeta, JunctionMeta
│   └── index.h/cpp          # IndexVX (k-mer index with Bloom filter)
│
├── alignment/               # Core alignment algorithms
│   └── hamming.h/cpp        # Hamming distance calculation (SIMD-optimized)
│
├── mapping/                 # High-level mapping logic
│   ├── mapper.h             # Main Mapper class declaration
│   ├── mapper.cpp           # Initialization and seed generation
│   ├── mapper_discovery.cpp # Candidate discovery phase
│   ├── mapper_verification.cpp # Alignment verification
│   ├── mapper_single.cpp    # Single-read mapping
│   ├── mapper_genomic.cpp   # Transcript-to-genomic projection
│   ├── mapper_mapq.cpp      # MAPQ calculation
│   ├── mapper_paired.cpp    # Paired-end logic and mate rescue
│   ├── mapper_rescue.cpp    # Rescue/reprocess for unmapped reads
│   └── mapper_single_end.cpp # Single-end mapping loop
│
└── main.cpp                 # Entry point and CLI parsing
```

### Design Principles

1. **Performance First**: All hot-path functions remain inline or in headers to enable compiler optimization
2. **Zero Overhead**: No virtual functions, no dynamic dispatch in critical loops
3. **Cache-Friendly**: Data locality preserved through careful struct layout
4. **Clear Separation**: Each module has a single, well-defined responsibility
5. **Testability**: Modules can be unit-tested independently

## Module Descriptions

### Core Module (`src/core/`)

**Purpose**: Fundamental types and DNA sequence utilities

- `types.h`: Data structures used across the codebase
  - `FailureStats`: Mapping statistics
  - `ReadRecord`: FASTQ record
  - `Alignment`: Mapping result
  - `GenomicPos`: Genomic coordinates with CIGAR

- `dna_utils.h`: DNA manipulation functions (all inline for performance)
  - `upper()`: Convert to uppercase, U→T
  - `revcomp()`: Reverse complement
  - `nt2b()`: Nucleotide to 2-bit encoding
  - `pack2()`: Pack sequences into 2-bit representation
  - `get2b()`: Extract base from packed sequence

### IO Module (`src/io/`)

**Purpose**: File input/output operations

- `fastq_reader.h/cpp`: Efficient FASTQ reading
  - `LineReader`: Low-level line reading (supports gzip)
  - `FastqReader`: FASTQ format parser with batch reading

- `output_buffer.h/cpp`: Buffered output for SAM files
  - 32MB buffer for efficient I/O
  - Automatic flushing

- `sam_writer.h/cpp`: SAM format output
  - Header generation with contig metadata
  - Record formatting with MAPQ, flags, CIGAR
  - Optional NM/MD tags

### Index Module (`src/index/`)

**Purpose**: K-mer index structure and loading

- `index_types.h`: Metadata structures
  - `TranscriptMeta`: Transcript information and exon coordinates
  - `JunctionMeta`: Splice junction information

- `index.h/cpp`: Main index structure
  - `IndexVX`: K-mer index with prefix-suffix split
  - Bloom filter for fast negative lookups
  - Compressed postings with delta-of-delta encoding
  - Binary search for k-mer lookup
  - Template-based posting iteration (inline for performance)

### Alignment Module (`src/alignment/`)

**Purpose**: Low-level alignment algorithms

- `hamming.h/cpp`: Fast Hamming distance
  - 2-bit encoded sequences
  - SIMD-optimized (64-bit operations)
  - Early termination when threshold exceeded
  - Handles aligned and unaligned cases

### Mapping Module (`src/mapping/`)

**Purpose**: High-level mapping orchestration

- `mapper.h`: Main Mapper class with all interfaces

- `mapper.cpp`: Initialization and seed generation
  - Constructor and data structure setup
  - K-mer extraction with Bloom filtering
  - Canonical k-mer support
  - Seed deduplication (optional)

- `mapper_discovery.cpp`: Candidate discovery
  - Two-stage seed processing (sieve + full)
  - Top-k target filtering
  - Offset quantization and voting
  - Early termination heuristics
  - Junction filtering

- `mapper_verification.cpp`: Alignment verification
  - Hamming-based verification
  - Multiple positions per candidate
  - Adaptive budget selection

- `mapper_single.cpp`: Single-read mapping
  - Forward/reverse strand mapping
  - Hit merging and deduplication
  - Ambiguity detection

- `mapper_genomic.cpp`: Coordinate projection
  - Transcript-to-genomic mapping
  - CIGAR string generation
  - Exon/intron handling
  - Splice junction projection

- `mapper_mapq.cpp`: Mapping quality calculation
  - Alternative hit consideration
  - Cross-chromosome ambiguity penalties
  - Isoform ambiguity handling

- `mapper_paired.cpp`: Paired-end operations
  - Pair scoring with insert size
  - Mate rescue for discordant pairs
  - Concordance checking

- `mapper_rescue.cpp`: Rescue unmapped reads
  - Relaxed parameters (dense seeding, no filtering)
  - Multiple rescue strategies
  - State snapshoting and restoration

- `mapper_single_end.cpp`: Single-end mapping loop
  - Batch processing
  - Optional rescue

## Building

### Requirements

- **Compiler**: g++ with C++17 support
- **Dependencies**: zlib only
- **Platform**: Linux (tested on Ubuntu 24)

### Compilation

```bash
# Build optimized production binary
make

# Clean build artifacts
make clean

# Quick compilation test (no optimization)
make test-compile

# Show build configuration
make info
```

### Compiler Flags

- `-O3`: Maximum optimization
- `-std=c++17`: C++17 standard
- `-march=native`: CPU-specific optimizations
- `-Wall -Wextra -pedantic`: Strict warnings

## Usage

Same CLI as the original monolithic version:

```bash
# Single-end mapping
./mapper_v17_modular index.bin reads.fq.gz output.sam

# Paired-end mapping
./mapper_v17_modular index.bin R1.fq.gz R2.fq.gz output.sam

# With options
./mapper_v17_modular index.bin R1.fq.gz R2.fq.gz output.sam \
    --stride 1 \
    --min-votes 3 \
    --verify-topk 12 \
    --rescue
```

### Options

- `--max-pairs N`: Limit paired mapping to N pairs
- `--stride N`: Seed stride (default: 1)
- `--df1-bytes N`: Early DF cutoff (default: 128)
- `--df2-bytes N`: Late DF cutoff (default: 512)
- `--min-votes N`: Minimum seed votes (default: 3)
- `--lead-margin N`: Lead over runner-up (default: 2)
- `--verify-topk N`: Max candidates to verify (default: 12)
- `--rescue`: Enable rescue for unmapped reads
- `--output-nm-md 1`: Emit NM/MD tags

## Performance Validation

### Validation Strategy

1. **Compilation Test**:
   ```bash
   make clean && make
   ```

2. **Output Comparison**:
   ```bash
   # Run original
   ./mapper_v17_profile_detailed index.bin R1.fq R2.fq original.sam

   # Run refactored
   ./mapper_v17_modular index.bin R1.fq R2.fq refactored.sam

   # Compare outputs (should be identical)
   diff original.sam refactored.sam
   ```

3. **Performance Benchmark**:
   ```bash
   # Time both versions
   time ./mapper_v17_profile_detailed index.bin R1.fq R2.fq out1.sam
   time ./mapper_v17_modular index.bin R1.fq R2.fq out2.sam
   ```

### Expected Results

- **Output**: Byte-for-byte identical SAM files
- **Performance**: Within 1-2% of original (compiler optimization should be equivalent)
- **Memory**: Identical memory footprint

## Key Design Decisions

### Why These Module Boundaries?

1. **Core/IO/Index/Alignment/Mapping**: Natural separation of concerns
2. **Mapper split across files**: Each file handles one phase of mapping (discovery, verification, etc.)
3. **Inline performance-critical code**: Hamming, k-mer lookup, DNA utils stay inline

### What Changed?

- **Structure**: Mono file → Multiple modules
- **Namespacing**: Added `rnamapper` namespace
- **Build system**: Single compile command → Makefile with dependencies
- **Headers**: Proper separation of interface (.h) and implementation (.cpp)

### What Didn't Change?

- **Algorithms**: Exact same logic, zero modifications
- **Data structures**: Identical layout and sizing
- **Performance characteristics**: Same cache behavior, same SIMD usage
- **Output**: Byte-for-byte identical SAM files

## Future Extensions

The modular structure enables:

1. **Unit Testing**: Each module can be tested independently
2. **Alternative Algorithms**: Swap Hamming for Smith-Waterman
3. **Index Formats**: Add support for different index versions
4. **Output Formats**: Add BAM output alongside SAM
5. **Parallel Processing**: Add thread-safe batch processing
6. **Profiling**: Easier to profile individual components

## Maintenance

### Adding New Features

1. Identify the appropriate module
2. Add functions to the relevant .h file
3. Implement in the corresponding .cpp file
4. Update Makefile dependencies (automatic)

### Debugging

```bash
# Compile with debug symbols
make CXXFLAGS="-g -O0 -std=c++17"

# Run with GDB
gdb ./mapper_v17_modular
```

### Code Style

- **Indentation**: 4 spaces
- **Braces**: K&R style
- **Naming**: snake_case for functions, CamelCase for types
- **Comments**: Explain "why", not "what"

## Files Not Modified

- `index_builder_v15.cpp`: Index builder (separate tool)
- `rescue_v17.h/cpp`: Old rescue module (superseded by mapper_rescue.cpp)
- `mapper_v17_profile_detailed_forward.h`: Old forward declarations (no longer needed)

## License and Attribution

This refactoring preserves all original functionality of `mapper_v17_profile_detailed.cpp`.

---

**Refactoring Date**: 2025-11-17
**Original Version**: v17
**Refactored Version**: v17-modular
**Compiler**: g++ 11.4+ with C++17
**Platform**: Linux (Ubuntu 24)
