include { MINIMAP2_ALIGN    as MINIMAP2_ALIGN_EXON6 } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN    as MINIMAP2_ALIGN_EXON7 } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_COVERAGE as BAM_COVERAGE_EXON6   } from '../../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_COVERAGE as BAM_COVERAGE_EXON7   } from '../../../modules/nf-core/samtools/coverage/main'

workflow MINIMAP2_ALIGN_READS {

    take:
    ch_samplesheet  // channel: [ val(meta), [ fastq ] ]
    ch_exon6fai     // channel: [ val(meta), [ fai ] ]
    ch_exon6fasta   // channel: [ val(meta), [ fasta ] ]
    ch_exon7fai     // channel: [ val(meta), [ fai ] ]
    ch_exon7fasta   // channel: [ val(meta), [ fasta ] ]

    main:

    ch_versions = Channel.empty()

    // 
    // MODULE: Minimap2/align
    // 
    MINIMAP2_ALIGN_EXON6 (
        ch_samplesheet,
        ch_exon6fasta,
        bam_format="bam",
        bam_index_extension="bai",
        cigar_paf_format=false,
        cigar_bam=false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_EXON6.out.versions.first())

    MINIMAP2_ALIGN_EXON7 (
        ch_samplesheet,
        ch_exon7fasta,
        bam_format="bam",
        bam_index_extension="bai",
        cigar_paf_format=false,
        cigar_bam=false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_EXON7.out.versions.first())

    //
    // MODULE: Samtools/coverage 
    // 
    BAM_COVERAGE_EXON6 (
        MINIMAP2_ALIGN_EXON6.out.bam
            .join(MINIMAP2_ALIGN_EXON6.out.index),
        ch_exon6fasta,
        ch_exon6fai
    )
    ch_versions = ch_versions.mix(BAM_COVERAGE_EXON6.out.versions.first())

    BAM_COVERAGE_EXON7 (
        MINIMAP2_ALIGN_EXON7.out.bam
            .join(MINIMAP2_ALIGN_EXON7.out.index),
        ch_exon7fasta,
        ch_exon7fai
    )
    ch_versions = ch_versions.mix(BAM_COVERAGE_EXON7.out.versions.first())

    
    emit:
    exon6bam      = MINIMAP2_ALIGN_EXON6.out.bam           // channel: [ val(meta), [ bam ] ]
    exon6bai      = MINIMAP2_ALIGN_EXON6.out.index          // channel: [ val(meta), [ bai ] ]
    exon7bam      = MINIMAP2_ALIGN_EXON7.out.bam           // channel: [ val(meta), [ bam ] ]
    exon7bai      = MINIMAP2_ALIGN_EXON7.out.index          // channel: [ val(meta), [ bai ] ]

    versions      = ch_versions                     // channel: [ versions.yml ]
}

