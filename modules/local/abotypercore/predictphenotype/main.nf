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
    // abotyper-core 0.4.1 is on PyPI; its Bioconda recipe is in review. Until the Seqera Containers image is built
    // from that recipe, docker/podman use a local image of the same release (assets/dev/Dockerfile.abotyper-core).
    container "abotyper.local/abotyper-core:0.4.1"

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
