/*
  SUBWORKFLOW: PREDICTABOPHENOTYPE
*/

//  Description: Predicts ABO phenotype by combining variant frequency data
//   with BAM coverage, extracting SNPs, and mapping them to phenotypes.
//   Every per-sample file the aggregation reads is staged into the task;
//   nothing is read back from params.outdir.

include { ABO_GETABOSNPS } from '../../../modules/local/abo/getabosnps/main'
include { ABO_SNPS2PHENO } from '../../../modules/local/abo/snps2pheno/main'

workflow PREDICTABOPHENOTYPE {
    take:
    ch_variants_freq   // channel: [ val(meta), [ freq ] ]
    ch_bam_coverage    // channel: [ val(meta), [ cov ] ]
    ch_haplotypes      // channel: [ val(meta), path(*.Haplotypes.tsv) ]          - from ABO_CLAIR2HAPLOTYPES
    ch_evidence_status // channel: [ val(meta), path(*.evidence_status.tsv) ]     - from ABO_EVIDENCE
    ch_concordance     // channel: [ val(meta), path(*.concordance.tsv) ]         - from ABO_EVIDENCE
    ch_bam_stats       // channel: [ val(meta), path(*.stats) ]                   - from SAMTOOLS_STATS
    ch_bam_flagstat    // channel: [ val(meta), path(*.flagstat) ]                - from SAMTOOLS_FLAGSTAT
    ch_panel           // channel: path(variant panel) - params.abo_panel, staged into each task; [] for the bundled panel

    main:

    // JOIN: Variant frequency with BAM coverage
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)

    /*
    MODULE: ABO_GETABOSNPS
    */
    ABO_GETABOSNPS(
        ch_combined_input,
        ch_panel
    )

    // Every per-sample file aggregate_abo_reports.py reads, staged into the task
    ch_sample_ids = ABO_GETABOSNPS.out.phenotype
        .map { meta, _phenotype -> meta.id }
        .collect()

    ch_sample_files = ABO_GETABOSNPS.out.phenotype
        .mix(ch_haplotypes, ch_evidence_status, ch_concordance, ch_bam_coverage, ch_bam_stats, ch_bam_flagstat)
        .map { _meta, sample_file -> sample_file }
        .collect()

    // Sample table from the samplesheet: id, control type, and the sample name and barcode to report
    // (meta.sample_name / meta.barcode, set when the samplesheet is read). One row per sample, so it is never empty.
    ch_samples = ABO_GETABOSNPS.out.phenotype
        .map { meta, _phenotype -> "${meta.id}\t${meta.control ?: ''}\t${meta.sample_name ?: ''}\t${meta.barcode ?: ''}\n" }
        .collectFile(name: 'samples.tsv', sort: true)

    /*
    MODULE: ABO_SNPS2PHENO
    */
    ABO_SNPS2PHENO(
        ch_sample_ids,
        ch_sample_files,
        ch_samples,
        ch_panel
    )

    emit:
    abo_results      = ABO_SNPS2PHENO.out.txt
    evidence_summary = ABO_SNPS2PHENO.out.mqc
}
