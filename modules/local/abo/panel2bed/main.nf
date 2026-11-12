/*
  MODULE: ABO_PANEL2BED
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
process ABO_PANEL2BED {
    tag "${panel_file.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/abotyper-core:0.1.0--pyhdfd78af_0'
        : 'biocontainers/abotyper-core:0.1.0--pyhdfd78af_0'}"

    input:
    path panel_file

    output:
    path "*.bed", emit: bed
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'abo_variant_panel'
    """
    panel_to_bed.py \\
        --panel ${panel_file} \\
        --output ${prefix}.bed \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: 'abo_variant_panel'
    """
    touch ${prefix}.bed
    """
}
