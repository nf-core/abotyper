/*
  MODULE: ABO_SNPS2PHENO
  Description:
    Runs aggregate_abo_reports.py over every sample's combined/ directory
    (*.ABOPhenotype.txt from ABO_GETABOSNPS, *.Haplotypes.tsv from
    ABO_CLAIR2HAPLOTYPES), scoring primary and named-marker calls per the
    panel's thresholds, resolving Phenotype/Genotype/ExtendedGenotype, and
    writing the pipeline's final output: ABO_result.txt/.xlsx and
    final_export.csv. Runs once per pipeline run, not per sample.
*/
process ABO_SNPS2PHENO {
    tag "COMPILING ABO RESULTS"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pandas_python_pyyaml_xlsxwriter:e206f33f53f3639b'
        : 'community.wave.seqera.io/library/pandas_python_pyyaml_xlsxwriter:4325a71f45ad8e44'}"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path samples_dir
    path per_sample_processing

    output:
    path "ABO_result.txt", emit: txt
    path "ABO_result.xlsx", emit: xls
    path "ABO_results.log", emit: log
    path "final_export.csv", emit: csv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('numpy'), eval('python -c "import numpy; print(numpy.__version__)"'), emit: versions_numpy, topic: versions
    tuple val("${task.process}"), val('xlsxwriter'), eval('python -c "import xlsxwriter; print(xlsxwriter.__version__)"'), emit: versions_xlsxwriter, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    aggregate_abo_reports.py \\
        ${per_sample_processing} ${args} 2>&1 | tee ABO_results.log
    """

    stub:
    """
    touch final_export.csv
    touch ABO_results.log
    touch ABO_result.xlsx
    touch ABO_result.txt
    """
}
