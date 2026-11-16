/*
  SUBWORKFLOW: PREDICTABOPHENOTYPE
*/

//  Description: Predicts ABO phenotype by combining variant frequency data
//   with BAM coverage, extracting SNPs, and mapping them to phenotypes.
//   Every per-sample file the aggregation reads is staged into the task;
//   nothing is read back from params.outdir.

include { ABOTYPERCORE_PREDICTPHENOTYPE } from '../../../modules/local/abotypercore/predictphenotype/main'
include { ABOTYPERCORE_AGGREGATE        } from '../../../modules/local/abotypercore/aggregate/main'

workflow PREDICTABOPHENOTYPE {
    take:
    ch_variants_freq   // channel: [ val(meta), [ freq ] ]
    ch_bam_coverage    // channel: [ val(meta), [ cov ] ]
    ch_haplotypes      // channel: [ val(meta), path(*.Haplotypes.tsv) ]          - from ABOTYPERCORE_CLAIR2HAPLOTYPES
    ch_evidence_status // channel: [ val(meta), path(*.evidence_status.tsv) ]     - from ABOTYPERCORE_PILEUP2METRICS
    ch_concordance     // channel: [ val(meta), path(*.concordance.tsv) ]         - from ABOTYPERCORE_PILEUP2METRICS
    ch_bam_stats       // channel: [ val(meta), path(*.stats) ]                   - from SAMTOOLS_STATS
    ch_bam_flagstat    // channel: [ val(meta), path(*.flagstat) ]                - from SAMTOOLS_FLAGSTAT
    ch_panel           // channel: path(variant panel) - params.abo_panel, staged into each task; [] for the bundled panel

    main:

    // JOIN: Variant frequency with BAM coverage
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)

    /*
    MODULE: ABOTYPERCORE_PREDICTPHENOTYPE
    */
    ABOTYPERCORE_PREDICTPHENOTYPE(
        ch_combined_input,
        ch_panel
    )

    // Each sample's files, joined by sample (a sample without Clair3 haplotypes gets null there)
    ch_per_sample = ABOTYPERCORE_PREDICTPHENOTYPE.out.phenotype
        .join(ch_haplotypes, remainder: true)
        .join(ch_evidence_status, remainder: true)
        .join(ch_concordance, remainder: true)
        .join(ch_bam_coverage, remainder: true)
        .join(ch_bam_stats, remainder: true)
        .join(ch_bam_flagstat, remainder: true)
        .filter { row -> row[1] != null }   // samples that reached phenotype prediction

    // Every file abotyper-aggregate reads, staged into the task under sample_files/
    ch_sample_files = ch_per_sample
        .flatMap { row -> row[1..-1].findAll { sample_file -> sample_file != null } }
        .collect()

    // Sample table: id, the samplesheet's control type, sample name and barcode (meta.sample_name / meta.barcode,
    // set when the samplesheet is read), and each sample's staged files. One row per sample, so it is never empty.
    def file_columns = ['phenotype', 'haplotypes', 'evidence_status', 'concordance', 'coverage', 'stats', 'flagstat']
    def header = (['id', 'control', 'sample', 'barcode'] + file_columns).join('\t') + '\n'
    ch_samples = ch_per_sample
        .map { row ->
            def meta  = row[0]
            def files = row[1..-1].collect { sample_file -> sample_file ? "sample_files/${sample_file.name}" : '' }
            ([meta.id, meta.control ?: '', meta.sample_name ?: '', meta.barcode ?: ''] + files).join('\t') + '\n'
        }
        .collectFile(name: 'samples.tsv', seed: header, sort: true)

    /*
    MODULE: ABOTYPERCORE_AGGREGATE
    */
    ABOTYPERCORE_AGGREGATE(
        ch_sample_files,
        ch_samples,
        ch_panel
    )

    emit:
    abo_results      = ABOTYPERCORE_AGGREGATE.out.txt
    evidence_summary = ABOTYPERCORE_AGGREGATE.out.mqc
}
