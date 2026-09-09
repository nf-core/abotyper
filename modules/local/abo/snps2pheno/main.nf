process ABO_SNPS2PHENO {
    tag "COMPILING ABO RESULTS"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14bab3c5ebaf4e44be2bbc6fa56108905b8aceb3ecbc85371516dcf1cdaafd3a/data'
        : 'community.wave.seqera.io/library/pandas_xlsxwriter:b231bcbbf11b41fd'}"

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
        per_sample_processing ${args} 2>&1 | tee ABO_results.log
    """

    stub:
    """
    touch final_export.csv
    touch ABO_results.log
    touch ABO_result.xlsx
    touch ABO_result.txt
    """
}
