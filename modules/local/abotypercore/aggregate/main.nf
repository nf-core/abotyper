/*
  MODULE: ABOTYPERCORE_AGGREGATE
  Description:
    Runs abotyper-aggregate on every sample's files: *.ABOPhenotype.txt from
    ABOTYPERCORE_PREDICTPHENOTYPE, *.Haplotypes.tsv from ABOTYPERCORE_CLAIR2HAPLOTYPES,
    the ABOTYPERCORE_PILEUP2METRICS tables and the samtools coverage/stats/flagstat files.
    The files are staged into the task. A samples table names each sample's files
    (staged paths) together with the samplesheet's control type, sample name and barcode.
    Scores primary and named-marker calls per the panel's thresholds, resolves
    Phenotype/Genotype/ExtendedGenotype, and writes the pipeline's final output:
    ABO_result.txt/.xlsx and final_export.csv. Runs once per pipeline run, not per sample.
*/
process ABOTYPERCORE_AGGREGATE {
    tag "COMPILING ABO RESULTS"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6335514e9473e3b78713d0155b0f7968ecefcf84be243029fcce1ed504264d68/data'
        : 'community.wave.seqera.io/library/numpy_openpyxl_pandas_pip_pruned:317ee119edadb518'}"

    input:
    path sample_files, stageAs: 'sample_files/*'
    path samples
    path panel_file

    output:
    path "ABO_result.txt", emit: txt
    path "ABO_result.xlsx", emit: xls
    path "ABO_results.log", emit: log
    path "final_export.csv", emit: csv
    path "ABO_run_qc.txt", emit: qc, optional: true
    path "ABO_evidence_summary_mqc.tsv", emit: mqc, optional: true
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def panel = panel_file ? "--panel ${panel_file}" : ''
    def args  = task.ext.args ?: ''
    """
    set -o pipefail
    abotyper-aggregate \\
        --samples ${samples} \\
        ${panel} \\
        --threads ${task.cpus} \\
        ${args} 2>&1 | tee ABO_results.log
    """

    stub:
    """
    touch ABO_result.txt
    touch final_export.csv
    touch ABO_results.log
    touch ABO_result.xlsx
    touch ABO_run_qc.txt
    touch ABO_evidence_summary_mqc.tsv
    """
}
