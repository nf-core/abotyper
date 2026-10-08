/*
  MODULE: ABO_CLAIR2HAPLOTYPES
  Description:
    Converts a Clair3 phased VCF into the same Haplotypes.tsv format
    pysam_haploscan.py (HAPLOSCAN) produces, so PREDICTABOPHENOTYPE's
    existing PhaseConfidence scoring works against Clair3 phasing.

    Clair3's phasing is block-level (GT + PS tag), not per-read -- rows
    are synthesized in proportion to each phase set's own allele depths,
    not literal individual reads. See bin/clair2haplotypes.py docstring.
*/
process ABO_CLAIR2HAPLOTYPES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyyaml:10a1a757f29eaf02'
        : 'community.wave.seqera.io/library/python_pyyaml:a1e09ed0f4856f89'}"

    input:
    tuple val(meta), path(phased_vcf)
    path panel_file

    output:
    tuple val(meta), path("*.Haplotypes.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pyyaml'), eval('python3 -c "import yaml; print(yaml.__version__)"'), emit: versions_pyyaml, topic: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    clair2haplotypes.py \\
        -i ${phased_vcf} \\
        -o ${prefix}.Haplotypes.tsv \\
        --panel ${panel_file} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Haplotypes.tsv
    """
}
