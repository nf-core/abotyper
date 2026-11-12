/*
  MODULE: ABO_CLAIR2HAPLOTYPES
  Description:
    Converts a Clair3 phased VCF into the same Haplotypes.tsv format
    pysam_haploscan.py (HAPLOSCAN) produces, so PREDICTABOPHENOTYPE's
    existing PhaseConfidence scoring works against Clair3 phasing.

    Clair3's phasing is block-level (GT + PS tag), not per-read -- rows
    are synthesized in proportion to each phase set's own allele depths,
    not literal individual reads. See abotyper-core (clair2haplotypes).
*/
process ABO_CLAIR2HAPLOTYPES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b7f81ff95f37792e4c40a969be296b6e2e1e99000c0f288c82316eaa11e5eda/data'
        : 'community.wave.seqera.io/library/numpy_openpyxl_pandas_pip_pruned:81dee64f9dbf35a2'}"

    input:
    tuple val(meta), path(phased_vcf)
    path panel_file

    output:
    tuple val(meta), path("*.Haplotypes.tsv"), emit: tsv
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    script:
    def panel  = panel_file ? "--panel ${panel_file}" : ''
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clair2haplotypes.py \\
        -i ${phased_vcf} \\
        -o ${prefix}.Haplotypes.tsv \\
        ${panel} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Haplotypes.tsv
    """
}
