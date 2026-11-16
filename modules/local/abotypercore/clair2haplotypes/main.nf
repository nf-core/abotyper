/*
  MODULE: ABOTYPERCORE_CLAIR2HAPLOTYPES
  Description:
    Converts a Clair3 phased VCF into the same Haplotypes.tsv format
    pysam_haploscan.py (HAPLOSCAN) produces, so PREDICTABOPHENOTYPE's
    existing PhaseConfidence scoring works against Clair3 phasing.

    Clair3's phasing is block-level (GT + PS tag), not per-read -- rows
    are synthesized in proportion to each phase set's own allele depths,
    not literal individual reads. See abotyper-core (clair2haplotypes).
*/
process ABOTYPERCORE_CLAIR2HAPLOTYPES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // abotyper-core 0.4.1 is on PyPI; its Bioconda recipe is in review. Until the Seqera Containers image is built
    // from that recipe, docker/podman use a local image of the same release (assets/dev/Dockerfile.abotyper-core).
    container "abotyper.local/abotyper-core:0.4.1"

    input:
    tuple val(meta), path(phased_vcf)
    path panel_file

    output:
    tuple val(meta), path("*.Haplotypes.tsv"), emit: tsv
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def panel  = panel_file ? "--panel ${panel_file}" : ''
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abotyper-clair2haplotypes \\
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
