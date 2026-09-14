/*
  MODULE: ABO_PANEL2BED
  Description:
    Generates a BED file of diagnostic positions from the ABO variant panel
    (abo_variant_panel.yaml), for use as Clair3's --bed_fn and for
    bcftools view -R filtering of its output.

    Indel/homopolymer/multi-offset panel rows are padded on both sides
    (see bin/panel_to_bed.py); plain SNP rows stay an exact 1bp interval.
    Runs once against the reference-wide panel file (no per-sample meta),
    matching DOWNLOAD_CLAIR3_MODEL's broadcast pattern in
    variant_calling_haploscan.
*/
process ABO_PANEL2BED {
    tag "${panel_file.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8d69e246c0a530fa88ba496bf8a62bd282e770fb4b68d2877ab777ce4943fed1/data'
        : 'community.wave.seqera.io/library/json5_pandas_python:e3184b0698afebbd'}"

    input:
    path panel_file

    output:
    path "*.bed", emit: bed
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pyyaml'), eval('python3 -c "import yaml; print(yaml.__version__)"'), emit: versions_pyyaml, topic: versions

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
