/*
 * Subworkflow: minimap_align_exons
 * Description: Aligns each sample against the single combined ABO reference
 * using Minimap2, then runs Samtools coverage, flagstat, and stats.
 * Single-reference mode: every sample is aligned exactly once, no
 * per-exon metadata matching required.
 */

include { MINIMAP2_ALIGN    } from '../../../modules/nf-core/minimap2/align'
include { SAMTOOLS_COVERAGE } from '../../../modules/nf-core/samtools/coverage'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats'

workflow MINIMAP2_ALIGN_READS {
    take:
    ch_samplesheet     // channel: [ val(meta), [ fastq ] ]
    ch_reference_fasta // channel: [ val(meta1), path(fasta) ]
    ch_reference_fai   // channel: [ val(meta1), path(fai) ]

    main:

    // Value channels: single reference reused for every sample
    def ch_ref_fasta = ch_reference_fasta.first()

    // Reference tuple [meta, fasta, fai] for coverage/stats,
    // pinned as a value channel with .first
    def ch_ref_bundle = ch_reference_fasta
        .combine(ch_reference_fai)
        .map { meta_fa, fa, meta_fai, fai -> [ [id: 'ABO_REF'], fa, fai ] }
        .first()

    //
    // MODULE: MINIMAP2_ALIGN
    //
    MINIMAP2_ALIGN(
        ch_samplesheet,
        ch_ref_fasta,
        "bam",
        "bai",
        false,
        false,
    )

    ch_bam_bai = MINIMAP2_ALIGN.out.bam.join(MINIMAP2_ALIGN.out.index, by: 0)

    //
    // MODULE: SAMTOOLS_COVERAGE
    //
    SAMTOOLS_COVERAGE(
        ch_bam_bai,
        ch_ref_bundle,
    )

    //
    // MODULE: SAMTOOLS_FLAGSTAT
    //
    SAMTOOLS_FLAGSTAT(
        ch_bam_bai
    )

    //
    // MODULE: SAMTOOLS_STATS
    //
    SAMTOOLS_STATS(
        ch_bam_bai,
        ch_ref_bundle,
    )

    emit:
    bam      = MINIMAP2_ALIGN.out.bam          // channel: [ val(meta), path(bam) ]
    bai      = MINIMAP2_ALIGN.out.index        // channel: [ val(meta), path(bai) ]
    coverage = SAMTOOLS_COVERAGE.out.coverage  // channel: [ val(meta), path(txt) ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat  // channel: [ val(meta), path(flagstat) ]
    stats    = SAMTOOLS_STATS.out.stats        // channel: [ val(meta), path(stats) ]
    fasta    = ch_reference_fasta              // channel: [ val(meta1), path(fasta) ]
    fai      = ch_reference_fai                // channel: [ val(meta1), path(fai) ]
}
