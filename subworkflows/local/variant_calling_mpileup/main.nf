include { SAMTOOLS_MPILEUP   } from '../../../modules/nf-core/samtools/mpileup/main'
include { MPILEUP_NUCL_FREQ  } from '../../../modules/local/mpileupstats/main'

workflow VARIANTS_QUANTIFICATION {

    take:
    ch_bam       // channel: [ val(meta), path(bam)]    - BAM file from minimap2 alignment
    ch_bai       // channel: [ val(meta), path(bai)]    - BAI file from minimap2 alignment
    ch_fasta     // channel: [ val(meta1), path(fasta)] - fasta file from params exon6 and params exon 7
    ch_fai       // channel: [ val(meta1), path(fai)]   - fasta index for the references
    ch_bed       // channel: [ val(meta1), path(bed)]   - bed file for the references

    main:

    ch_versions = Channel.empty()

    // Prepare mpileup input by matching bam with bed files based on exon metadata
    ch_mpileup_input = ch_bam
        .combine(ch_bed)
        .filter { bam_meta, bam, bed_meta, bed ->
            bam_meta.exon == bed_meta.exon &&
            bam_meta.ref_id == bed_meta.id.replace('.fasta', '')
        }
        .map { bam_meta, bam, bed_meta, bed ->
            [bam_meta, bam, bed]
        }

    // Match fasta files with samples based on exon and reference ID
    ch_fasta_matched = ch_mpileup_input
        .combine(ch_fasta)
        .filter { bam_meta, bam, bed, fasta_meta, fasta ->
            bam_meta.exon == fasta_meta.exon &&
            bam_meta.ref_id == fasta_meta.id
        }
        .map { bam_meta, bam, bed, fasta_meta, fasta ->
            [bam_meta, bam, bed, fasta_meta, fasta]
        }

    /*
    MODULE: SAMTOOLS_MPILEUP
    */
    SAMTOOLS_MPILEUP (
        ch_fasta_matched.map { bam_meta, bam, bed, fasta_meta, fasta -> [bam_meta, bam, bed] },
        ch_fasta_matched.map { bam_meta, bam, bed, fasta_meta, fasta -> [fasta_meta, fasta] }
    )

    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

    // Match mpileup output with fasta for nucleotide frequency analysis
    ch_nucl_freq_input = SAMTOOLS_MPILEUP.out.mpileup
        .combine(ch_fasta)
        .filter { mpileup_meta, mpileup, fasta_meta, fasta ->
            mpileup_meta.exon == fasta_meta.exon
        }

    /*
    MODULE: MPILEUP_NUCL_FREQ
    */
    MPILEUP_NUCL_FREQ (
        ch_nucl_freq_input.map { mpileup_meta, mpileup, fasta_meta, fasta -> [mpileup_meta, mpileup] },
        ch_nucl_freq_input.map { mpileup_meta, mpileup, fasta_meta, fasta -> [fasta_meta, fasta] }
    )

    ch_versions = ch_versions.mix(MPILEUP_NUCL_FREQ.out.versions.first())

    emit:
    metrics       = MPILEUP_NUCL_FREQ.out.tsv              // channel: [ val(meta), path(tsv) ]
    versions      = ch_versions                            // channel: [ versions.yml ]
}
