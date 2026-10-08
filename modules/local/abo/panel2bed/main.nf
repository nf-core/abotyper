/*
  MODULE: ABO_PANEL2BED
  Description:
    Generates a BED file of diagnostic positions from the ABO variant panel
    (abo_variant_panel.yaml), for the VARIANTS_QC subworkflow's
    bcftools view -R filtering. Clair3 itself is deliberately NOT restricted
    to this BED via --bed_fn -- that would mean editing the vendored
    modules/nf-core/clair3 module, which is off-limits; every result this
    pipeline has been validated against ran Clair3 unrestricted.

    Indel/homopolymer/multi-offset panel rows are padded on both sides
    (see bin/panel_to_bed.py); plain SNP rows stay an exact 1bp interval.
    Runs once against the reference-wide panel file (no per-sample meta),
    matching DOWNLOAD_CLAIR3_MODEL's broadcast pattern in
    variants_quantification.
*/
process ABO_PANEL2BED {
    tag "${panel_file.baseName}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyyaml:10a1a757f29eaf02'
        : 'community.wave.seqera.io/library/python_pyyaml:a1e09ed0f4856f89'}"

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
