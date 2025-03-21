include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_EXON6   } from '../../../modules/nf-core/samtools/mpileup/main'
include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_EXON7   } from '../../../modules/nf-core/samtools/mpileup/main'
include { MAKEINDEX                                    } from '../../../modules/local/makeindex/main'
include { MPILEUP_NUCL_FREQ  as MPILEUP_EXON6_NUCL_FREQ  } from '../../../modules/local/mpileupstats/main'
include { MPILEUP_NUCL_FREQ  as MPILEUP_EXON7_NUCL_FREQ  } from '../../../modules/local/mpileupstats/main'

workflow VARIANTS_QUANTIFICATION {

    take:
    exon6bai     // channel: [ val(meta), [ bai ] ]
    exon6bam     // channel: [ val(meta), [ bam ] ]
    exon7bai     // channel: [ val(meta), [ bai ] ]
    exon7bam     // channel: [ val(meta), [ bam ] ]
    exon6fai     // channel: [ val(meta), [ fai ] ]
    exon6fasta   // channel: [ val(meta), [ fasta ]]
    exon7fai     // channel: [ val(meta), [ fai ] ]
    exon7fasta   // channel: [ val(meta), [ fasta ]]


    main:

    ch_versions = Channel.empty()

    //
    // MODULE: makeindex
    //
    MAKEINDEX (
        exon6fai,
        exon7fai
    )

    /*
    Some sanity check to ensure paths and metadata are as expected by the processes below. 
    SAMtools is the main offender. Modules parameterization keeps fluctuating between samtools/subtool 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // @TODO NOTE TO SELF:  UNCOMMENT TO EXECUTE FOR TESTING
    // exon6bam.view{ "exon6bam: $it" }
    // MAKEINDEX.out.exon6bed.view{ "MAKEINDEX.out.exon6bed: $it" }
    // exon6fasta.view{ "exon6fasta: $it" }

    // exon6bam.combine(MAKEINDEX.out.exon6bed)
    //     .map { bam_meta, bam, bed_meta, bed -> [bam_meta, bam, bed] }
    //     .view { "MPILEUP_EXON6 input: $it" }

    // exon7bam.view{ "exon7bam: $it" }
    // MAKEINDEX.out.exon7bed.view{ "MAKEINDEX.out.exon7bed: $it" }
    // exon7fasta.view{ "exon7fasta: $it" }

    // exon7bam.combine(MAKEINDEX.out.exon7bed)
    //     .map { bam_meta, bam, bed_meta, bed -> [bam_meta, bam, bed] }
    //     .view { "MPILEUP_EXON7 input: $it" }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    //
    // MODULE: samtools/mpileup
    //
    SAMTOOLS_MPILEUP_EXON6 (
        exon6bam.combine(MAKEINDEX.out.exon6bed.map{ meta, bed -> bed }),
        exon6fasta.map { meta, fasta -> fasta }
    )
    
    SAMTOOLS_MPILEUP_EXON7 ( 
        exon7bam.combine (MAKEINDEX.out.exon7bed.map{meta, bed -> bed}),
        exon7fasta.map { meta, fasta -> fasta }
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP_EXON6.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP_EXON7.out.versions.first())

    // 
    // MODULE: mpileupmetrics
    // 
    MPILEUP_EXON6_NUCL_FREQ (
        SAMTOOLS_MPILEUP_EXON6.out.mpileup,
        exon6fasta
    )

    MPILEUP_EXON7_NUCL_FREQ (
        SAMTOOLS_MPILEUP_EXON7.out.mpileup,
        exon7fasta
    )

    ch_versions = ch_versions.mix(MPILEUP_EXON6_NUCL_FREQ.out.versions.first())
    ch_versions = ch_versions.mix(MPILEUP_EXON7_NUCL_FREQ.out.versions.first())

    emit:
    exon6metrics      = MPILEUP_EXON6_NUCL_FREQ.out.tsv     // channel: [ val(meta), [ csv ] ]
    exon7metrics      = MPILEUP_EXON7_NUCL_FREQ.out.tsv     // channel: [ val(meta), [ csv ] ]
    versions          = ch_versions                     // channel: [ versions.yml ]
}

