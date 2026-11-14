/*
  SUBWORKFLOW: PREDICTABOPHENOTYPE
*/

//  Description: Predicts ABO phenotype by combining variant frequency data
//   with BAM coverage, extracting SNPs, and mapping them to phenotypes, then
//   aggregates every sample: the read-level call (ABOTYPERCORE_READTYPER, every sample)
//   with the per-position call (only the samples --clair3_mode sent to Clair3 + pileup).
//   Every per-sample file the aggregation reads is staged into the task;
//   nothing is read back from params.outdir.

include { ABOTYPERCORE_PREDICTPHENOTYPE } from '../../../modules/local/abotypercore/predictphenotype/main'
include { ABOTYPERCORE_AGGREGATE        } from '../../../modules/local/abotypercore/aggregate/main'

workflow PREDICTABOPHENOTYPE {
    take:
    ch_variants_freq   // channel: [ val(meta), [ freq ] ]                                - samples with per-position typing
    ch_bam_coverage    // channel: [ val(meta), [ cov ] ]                                 - every aligned sample
    ch_haplotypes      // channel: [ val(meta), path(*.Haplotypes.tsv) ]                  - from ABOTYPERCORE_CLAIR2HAPLOTYPES
    ch_evidence_status // channel: [ val(meta), path(*.evidence_status.tsv) ]             - from ABOTYPERCORE_PILEUP2METRICS
    ch_concordance     // channel: [ val(meta), path(*.concordance.tsv) ]                 - from ABOTYPERCORE_PILEUP2METRICS
    ch_bam_stats       // channel: [ val(meta), path(*.stats) ]                           - from SAMTOOLS_STATS
    ch_bam_flagstat    // channel: [ val(meta), path(*.flagstat) ]                        - from SAMTOOLS_FLAGSTAT
    ch_readtyping      // channel: [ val(meta), path(*.ReadTyping.tsv) ]                  - from ABOTYPERCORE_READTYPER
    ch_readtyping_haps // channel: [ val(meta), path(*.ReadTypingHaplotypes.tsv) ]        - from ABOTYPERCORE_READTYPER
    ch_panel           // channel: path(variant panel) - params.abo_panel, staged into each task; [] for the bundled panel

    main:

    // JOIN: Variant frequency with BAM coverage (samples sent to the per-position caller only)
    ch_combined_input = ch_variants_freq
        .join(ch_bam_coverage)

    /*
    MODULE: ABOTYPERCORE_PREDICTPHENOTYPE
    */
    ABOTYPERCORE_PREDICTPHENOTYPE(
        ch_combined_input,
        ch_panel
    )

    // Each sample's files, joined by sample. A file a sample does not have is null: the Clair3 haplotypes when
    // Clair3 failed, and every per-position file when --clair3_mode did not send the sample to the per-position
    // caller. The columns follow the join order; they are named in the table header, so the order does not matter.
    def file_columns = ['coverage', 'phenotype', 'haplotypes', 'evidence_status', 'concordance', 'stats', 'flagstat',
                        'readtyping', 'readtyping_haplotypes']
    ch_per_sample = ch_bam_coverage
        .join(ABOTYPERCORE_PREDICTPHENOTYPE.out.phenotype, remainder: true)
        .join(ch_haplotypes, remainder: true)
        .join(ch_evidence_status, remainder: true)
        .join(ch_concordance, remainder: true)
        .join(ch_bam_stats, remainder: true)
        .join(ch_bam_flagstat, remainder: true)
        .join(ch_readtyping, remainder: true)
        .join(ch_readtyping_haps, remainder: true)
        .filter { row -> row[2] != null || row[8] != null }   // typed by at least one method

    // Every file abotyper-aggregate reads, staged into the task under sample_files/
    ch_sample_files = ch_per_sample
        .flatMap { row -> row[1..-1].findAll { sample_file -> sample_file != null } }
        .collect()

    // Sample table: id, the samplesheet's control type, sample name, barcode and expected genotype (positive
    // controls; meta.sample_name / meta.barcode / meta.expected, set when the samplesheet is read), and each
    // sample's staged files. One row per sample, so it is never empty.
    def header = (['id', 'control', 'sample', 'barcode', 'expected'] + file_columns).join('\t') + '\n'
    ch_samples = ch_per_sample
        .map { row ->
            def meta  = row[0]
            def files = row[1..-1].collect { sample_file -> sample_file ? "sample_files/${sample_file.name}" : '' }
            ([meta.id, meta.control ?: '', meta.sample_name ?: '', meta.barcode ?: '', meta.expected ?: ''] + files)
                .join('\t') + '\n'
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
    readtyping_mqc   = ABOTYPERCORE_AGGREGATE.out.readtyping_mqc
    novel_variants   = ABOTYPERCORE_AGGREGATE.out.novel_variants
}
