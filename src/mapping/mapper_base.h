// mapper_base.h - Abstract base class for dual-mode RNA-seq mapping
// Provides common interface for both Illumina (paired-end) and ONT (long-read) mappers
#pragma once

#include "../core/types.h"
#include "../index/index.h"
#include <vector>
#include <string>
#include <cstdint>

namespace rnamapper {

/**
 * MapperBase - Abstract base class for RNA-seq mappers
 *
 * This class defines the common interface that both IlluminaMapper and ONTMapper
 * must implement. It provides:
 * - Shared access to the reference index (IndexVX)
 * - Shared statistics tracking (FailureStats)
 * - Virtual methods for the three main mapping modes
 *
 * Design rationale:
 * - Illumina reads are short (50-300bp), paired-end, low error rate (~1%)
 * - ONT reads are long (1kb-100kb), single-molecule, higher error rate (~5-15%)
 * - Different mapping strategies are needed, but the interface remains the same
 */
class MapperBase {
protected:
    // Reference to the k-mer index (shared by all mapper implementations)
    const IndexVX& IX;

public:
    // Statistics tracking (shared by all mapper implementations) - public for SAMWriter access
    FailureStats& ST;
    /**
     * Constructor - Initialize mapper with index and statistics
     * @param ix Reference to the loaded k-mer index
     * @param st Reference to statistics structure for tracking failures/successes
     */
    MapperBase(const IndexVX& ix, FailureStats& st) : IX(ix), ST(st) {}

    /**
     * Virtual destructor - Required for proper cleanup of derived classes
     */
    virtual ~MapperBase() = default;

    // ==================== CORE MAPPING INTERFACE ====================

    /**
     * Map a single read sequence to the transcriptome
     *
     * This is the core single-read mapping function. Given a read name and sequence,
     * it returns all valid alignments found.
     *
     * @param qname Query name (read identifier from FASTQ)
     * @param seq Read sequence (DNA string)
     * @return Vector of alignments (may be empty if no valid alignments found)
     *
     * Implementation notes:
     * - Illumina: Uses k-mer seeding with tight mismatch thresholds
     * - ONT: May use minimizers or sparse seeding with relaxed thresholds
     */
    virtual std::vector<Alignment> map_read_single(
        const std::string& qname,
        const std::string& seq) = 0;

    /**
     * Map paired-end reads from two FASTQ files
     *
     * This function processes paired-end sequencing data, finding concordant
     * alignments for read pairs.
     *
     * @param r1_path Path to R1 FASTQ file (forward reads)
     * @param r2_path Path to R2 FASTQ file (reverse reads)
     * @param out_path Path to output SAM file
     * @param max_pairs Maximum number of pairs to process (0 = unlimited)
     * @param do_rescue Enable mate rescue for orphaned reads
     *
     * Implementation notes:
     * - Illumina: Primary use case, includes insert size filtering and mate rescue
     * - ONT: Not typically used (ONT is single-molecule), may throw or return early
     */
    virtual void map_paired(
        const std::string& r1_path,
        const std::string& r2_path,
        const std::string& out_path,
        uint64_t max_pairs = 0,
        bool do_rescue = false) = 0;

    /**
     * Map single-end reads from a FASTQ file
     *
     * This function processes single-end sequencing data, mapping each read
     * independently.
     *
     * @param reads_path Path to FASTQ file
     * @param out_path Path to output SAM file
     * @param max_reads Maximum number of reads to process (0 = unlimited)
     * @param do_rescue Enable reprocessing for failed reads (Illumina-specific)
     *
     * Implementation notes:
     * - Illumina: Used for single-end libraries, less common than paired-end
     * - ONT: Primary use case, processes long single-molecule reads
     */
    virtual void map_single_end(
        const std::string& reads_path,
        const std::string& out_path,
        uint64_t max_reads = 0,
        bool do_rescue = false) = 0;

    // ==================== OPTIONAL HOOKS ====================

    /**
     * Called at the start of each read/pair processing epoch
     *
     * This allows implementations to reset per-read state, update internal
     * data structures, or perform other housekeeping tasks.
     *
     * Default: No-op (can be overridden by derived classes if needed)
     */
    virtual void next_read_epoch() {}

    // ==================== ILLUMINA-SPECIFIC INTERFACE ====================
    // These methods are required by SAMWriter and other downstream components
    // They may not be relevant for all mapper types (e.g., ONT)

    /**
     * Convert transcript coordinates to genomic coordinates
     * @param tid Transcript ID
     * @param tx_pos Position in transcript
     * @param read_len Read length
     * @return Genomic position with CIGAR string
     */
    virtual GenomicPos transcript_to_genomic(uint32_t tid, int tx_pos, int read_len) const = 0;

    /**
     * Calculate mapping quality (MAPQ) for an alignment
     * @param best Best alignment (can be nullptr)
     * @param all_hits All candidate alignments
     * @param read_len Read length
     * @return MAPQ value (0-255)
     */
    virtual uint8_t calculate_mapq(const Alignment* best, const std::vector<Alignment>& all_hits, int read_len) const = 0;

    /**
     * Get minimum insert size for paired-end reads
     * @return Minimum insert size (Illumina-specific)
     */
    virtual int get_min_ins() const = 0;

    /**
     * Get maximum insert size for paired-end reads
     * @return Maximum insert size (Illumina-specific)
     */
    virtual int get_max_ins() const = 0;
};

} // namespace rnamapper
