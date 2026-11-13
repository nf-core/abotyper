/*
  MODULE: ABOTYPERCORE_PREDICTPHENOTYPE
  Description:
    Runs abotyper-predict-phenotype: reads the per-position metrics
    (AlignmentStatistics.tsv, from ABOTYPERCORE_CLAIR2METRICS) and BAM coverage for
    one sample, and writes a human-readable *.ABOPhenotype.txt report --
    one marker per diagnostic panel position, with its interpretation text.
    No -e/--exon is passed: the script auto-detects every exon present via
    panel position overlap and writes all sections into one combined report.
    This is the report ABOTYPERCORE_AGGREGATE later aggregates across samples.
*/
process ABOTYPERCORE_PREDICTPHENOTYPE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/6335514e9473e3b78713d0155b0f7968ecefcf84be243029fcce1ed504264d68/data'
        : 'community.wave.seqera.io/library/numpy_openpyxl_pandas_pip_pruned:317ee119edadb518'}"

    input:
    tuple val(meta), path(variants_freq), path(coverage)
    path panel_file

    output:
    tuple val(meta), path("*.ABOPhenotype.txt"), emit: phenotype
    tuple val(meta), path("*.log.txt"), emit: log
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def panel  = panel_file ? "--panel ${panel_file}" : ''
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    abotyper-predict-phenotype \\
        -i ${variants_freq} \\
        -o ${prefix}.ABOPhenotype.txt \\
        -c ${coverage} \\
        ${panel} \\
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
