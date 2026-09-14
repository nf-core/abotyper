/*
https://github.com/pjedge/longshot
Longshot is a variant calling tool for diploid genomes using long error prone reads such as
Pacific Biosciences (PacBio) SMRT and Oxford Nanopore Technologies (ONT).
It takes as input an aligned BAM/CRAM file and outputs a phased VCF file with variants and haplotype information.
It can also genotype and phase input VCF files. It can output
haplotype-separated BAM files that can be used for downstream analysis.
Currently, it only calls single nucleotide variants (SNVs),
but it can genotype indels if they are given in an input VCF.
*/

process LONGSHOT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longphase:1.7.3--hf5e1c6e_0':
        'biocontainers/longphase:1.7.3--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta1), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf"),      emit: vcf, optional: true
    path "versions.yml",                 emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    longshot \\
        $args \\
        --ref ${fasta} \\
        --out ${prefix}.longshot.vcf \\
        --bam $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longshot --version | awk -F " " '{print \$NF}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.longshot.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longshot --version | awk -F " " '{print \$NF}')
    END_VERSIONS
    """
}
