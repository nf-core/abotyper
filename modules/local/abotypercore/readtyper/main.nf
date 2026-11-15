/*
  MODULE: ABOTYPERCORE_READTYPER
  Description:
    Read-level ABO typing (abotyper-readtyper). Every read of the exon 2-7 amplicon is scored against every ISBT
    ABO allele, and the allele pair that best explains all reads is called, with its likelihood margin. Checks
    reported with the call: allele balance, reads that fit neither called allele (a third haplotype, possible
    contamination), chimeras of the two alleles and the minimum number of informative reads. The reads of each
    allele are then pooled into a haplotype consensus, which confirms rare alleles and places coding variants that
    are not in the ISBT table, and known sites where the haplotype differs from its allele, on the right allele.
    Needs only the BAM and the reference: no Clair3, no pileup.

    See abotyper-core (readtyper).
*/
process ABOTYPERCORE_READTYPER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // development image: abotyper-core 0.4.0.dev3 is not on Bioconda yet, so there is no Seqera Containers image.
    // Build it with assets/dev/Dockerfile.abotyper-core-dev (see docs/usage.md), or run with -profile conda.
    container "abotyper.local/abotyper-core:0.4.0.dev3"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.ReadTyping.tsv"),           emit: tsv
    tuple val(meta), path("*.ReadTypingHaplotypes.tsv"), emit: haplotypes
    tuple val("${task.process}"), val('abotyper-core'), eval('python3 -c "import abotyper_core; print(abotyper_core.__version__)"'), emit: versions_abotyper_core, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    abotyper-readtyper \\
        --bam ${bam} \\
        --ref ${fasta} \\
        --sample "${meta.id}" \\
        --out ${prefix}.ReadTyping.tsv \\
        --haplotypes-out ${prefix}.ReadTypingHaplotypes.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ReadTyping.tsv ${prefix}.ReadTypingHaplotypes.tsv
    """
}
