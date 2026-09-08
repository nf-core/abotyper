/*
  SUBWORKFLOW: VARIANTS_QUANTIFICATION (haploscan)

  Uses pysam_haploscan.py to compute per-position allele frequencies AND
  per-read haplotypes in a single pass over the aligned BAM — no pileup file.

  Single-reference mode: every sample is aligned against the same combined
  ABO reference, so no per-exon metadata matching is required — the FASTA
  and FAI are broadcast as value channels to every sample.
*/

include { HAPLOSCAN } from '../../../modules/local/haploscan/main'

workflow VARIANTS_QUANTIFICATION {
    take:
    ch_bam   // channel: [ val(meta),  path(bam)   ]
    ch_bai   // channel: [ val(meta),  path(bai)   ]
    ch_fasta // channel: [ val(meta1), path(fasta) ] - single combined reference
    ch_fai   // channel: [ val(meta1), path(fai)   ] - single combined reference

    main:

    // Join BAM and BAI (metadata matches exactly)
    ch_bam_bai = ch_bam.join(ch_bai, by: 0)

    // .combine always returns a queue channel. Adding .first() operator returns a reusable value channel.
    def ch_fasta_fai = ch_fasta
        .combine(ch_fai)
        .map { fasta_meta, fasta, fai_meta, fai -> [fasta_meta, fasta, fai] }
        .first()

    //
    // MODULE HAPLOSCAN
    //
    HAPLOSCAN(
        ch_bam_bai,
        ch_fasta_fai,
    )

    emit:
    metrics    = HAPLOSCAN.out.tsv        // channel: [ val(meta), path(*.AlignmentStatistics.tsv) ]
    haplotypes = HAPLOSCAN.out.haplotypes // channel: [ val(meta), path(*.Haplotypes.tsv) ]
}
