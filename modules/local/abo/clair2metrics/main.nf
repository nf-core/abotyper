/*
  MODULE: ABO_CLAIR2METRICS
  Description:
    Converts a Clair3 --gvcf output into the same AlignmentStatistics.tsv
    format pysam_haploscan.py (HAPLOSCAN) produces, so ABO_GETABOSNPS /
    predict_abo_phenotype.py need no changes to consume Clair3 calls.

    Only the panel's calibrated diagnostic positions are emitted. Works
    directly off Clair3's raw gVCF (no bcftools norm step required --
    validated against a real cohort; see the abotyper-core clair2metrics docstring).
    gVCF non-variant blocks supply real reference-depth coverage where no
    variant was called, which is what recovers indel signal (e.g.
    c.1061delC) that Clair3's plain merge_output silently drops.
*/
process ABO_CLAIR2METRICS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b7f81ff95f37792e4c40a969be296b6e2e1e99000c0f288c82316eaa11e5eda/data'
        : 'community.wave.seqera.io/library/numpy_openpyxl_pandas_pip_pruned:81dee64f9dbf35a2'}"

    input:
    tuple val(meta), path(gvcf)
    path panel_file

    output:
    tuple val(meta), path("*.AlignmentStatistics.tsv"), emit: tsv
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

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
