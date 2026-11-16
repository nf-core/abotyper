/*
  MODULE: ABOTYPERCORE_PANELTOBED
  Description:
    Generates a BED file of diagnostic positions from the ABO variant panel
    (abo_variant_panel.yaml), for the VARIANTS_QC subworkflow's
    bcftools view -R filtering. Clair3 is not restricted to this BED: it
    calls over the whole exon 2-7 region (params.abo_clair3_region_bed).

    Indel/homopolymer/multi-offset panel rows are padded on both sides
    (see abotyper-core, panel_to_bed); plain SNP rows stay an exact 1bp interval.
    Runs once against the reference-wide panel file (no per-sample meta),
    matching CLAIR3_DOWNLOADMODEL's broadcast pattern in
    variants_quantification.
*/
process ABOTYPERCORE_PANELTOBED {
    tag "${panel_file ? panel_file.baseName : 'abotyper-core panel'}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // abotyper-core 0.4.1 is on PyPI; its Bioconda recipe is in review. Until the Seqera Containers image is built
    // from that recipe, docker/podman use a local image of the same release (assets/dev/Dockerfile.abotyper-core).
    container "abotyper.local/abotyper-core:0.4.1"

    input:
    path panel_file

    output:
    path "*.bed", emit: bed
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def panel  = panel_file ? "--panel ${panel_file}" : ''
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'abo_variant_panel'
    """
    abotyper-panel-to-bed \\
        ${panel} \\
        --output ${prefix}.bed \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: 'abo_variant_panel'
    """
    touch ${prefix}.bed
    """
}
