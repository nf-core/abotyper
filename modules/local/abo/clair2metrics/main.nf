/*
  MODULE: ABO_CLAIR2METRICS
  Description:
    Converts a Clair3 --gvcf output into the same AlignmentStatistics.tsv
    format pysam_haploscan.py (HAPLOSCAN) produces, so ABO_GETABOSNPS /
    predict_abo_phenotype.py need no changes to consume Clair3 calls.

    Only the panel's calibrated diagnostic positions are emitted. Works
    directly off Clair3's raw gVCF (no bcftools norm step required --
    validated against a real cohort; see bin/clair2metrics.py docstring).
    gVCF non-variant blocks supply real reference-depth coverage where no
    variant was called, which is what recovers indel signal (e.g.
    c.1061delC) that Clair3's plain merge_output silently drops.
*/
process ABO_CLAIR2METRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyyaml:10a1a757f29eaf02'
        : 'community.wave.seqera.io/library/python_pyyaml:a1e09ed0f4856f89'}"

    input:
    tuple val(meta), path(gvcf)
    path panel_file

    output:
    tuple val(meta), path("*.AlignmentStatistics.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pyyaml'), eval('python3 -c "import yaml; print(yaml.__version__)"'), emit: versions_pyyaml, topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clair2metrics.py \\
        -i ${gvcf} \\
        -o ${prefix}.AlignmentStatistics.tsv \\
        --panel ${panel_file} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.AlignmentStatistics.tsv
    """
}
