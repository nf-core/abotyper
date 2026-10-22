/*
  MODULE: ABO_GETABOSNPS
  Description:
    Runs predict_abo_phenotype.py: reads the per-position metrics
    (AlignmentStatistics.tsv, from ABO_CLAIR2METRICS) and BAM coverage for
    one sample, and writes a human-readable *.ABOPhenotype.txt report --
    one marker per diagnostic panel position, with its interpretation text.
    No -e/--exon is passed: the script auto-detects every exon present via
    panel position overlap and writes all sections into one combined report.
    This is the report ABO_SNPS2PHENO later aggregates across samples.
*/
process ABO_GETABOSNPS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pandas_python_pyyaml:a1b47a7a3bb7f9c7'
        : 'community.wave.seqera.io/library/pandas_python_pyyaml:e40312b9d861eff3'}"

    input:
    tuple val(meta), path(variants_freq), path(coverage)

    output:
    tuple val(meta), path("*.ABOPhenotype.txt"), emit: phenotype
    tuple val(meta), path("*.log.txt"), emit: log
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    predict_abo_phenotype.py \\
        -i ${variants_freq} \\
        -o ${prefix}.ABOPhenotype.txt \\
        -c ${coverage} \\
        ${args} \\
        2>&1 | tee ${prefix}.log.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ABOPhenotype.txt
    touch ${prefix}.log.txt
    """
}
