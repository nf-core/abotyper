/*
  MODULE: ABO_SNPS2PHENO
  Description:
    Runs aggregate_abo_reports.py over every sample's files (*.ABOPhenotype.txt
    from ABO_GETABOSNPS, *.Haplotypes.tsv from ABO_CLAIR2HAPLOTYPES, the
    ABO_EVIDENCE tables and the samtools coverage/stats/flagstat files), which
    are staged into the task and laid out as per_sample_processing/<id>/combined/.
    Controls declared in the samplesheet arrive as a two-column table.
    Scores primary and named-marker calls per the panel's thresholds, resolves
    Phenotype/Genotype/ExtendedGenotype, and writes the pipeline's final output:
    ABO_result.txt/.xlsx and final_export.csv. Runs once per pipeline run, not per sample.
*/
process ABO_SNPS2PHENO {
    tag "COMPILING ABO RESULTS"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/abotyper-core:0.1.0--pyhdfd78af_0'
        : 'biocontainers/abotyper-core:0.1.0--pyhdfd78af_0'}"

    input:
    val sample_ids
    path sample_files, stageAs: 'sample_files/*'
    path controls
    path panel_file

    output:
    path "ABO_result.txt", emit: txt
    path "ABO_result.xlsx", emit: xls
    path "ABO_results.log", emit: log
    path "final_export.csv", emit: csv
    path "ABO_run_qc.txt", emit: qc, optional: true
    path "ABO_evidence_summary_mqc.tsv", emit: mqc, optional: true
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    script:
    def args = task.ext.args ?: ''
    def controls_arg = controls ? "--controls ${controls}" : ''
    // Build the per_sample_processing/<id>/combined/ layout the script reads, inside the task
    """
    for id in ${sample_ids.join(' ')}; do
        mkdir -p "per_sample_processing/\${id}/combined"
        for suffix in ABOPhenotype.txt Haplotypes.tsv evidence_status.tsv concordance.tsv coverage.txt stats flagstat; do
            if [ -e "sample_files/\${id}.\${suffix}" ]; then
                ln -s "../../../sample_files/\${id}.\${suffix}" "per_sample_processing/\${id}/combined/"
            fi
        done
    done

    set -o pipefail
    aggregate_abo_reports.py \\
        per_sample_processing \\
        --panel ${panel_file} \\
        ${controls_arg} \\
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
