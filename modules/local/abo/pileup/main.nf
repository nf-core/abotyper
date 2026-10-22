/*
  MODULE: ABO_PILEUP
  Description:
    samtools mpileup over the ABO panel BED (from ABO_PANEL2BED), straight from the
    aligned BAM. The second, independent line of evidence for every panel position,
    compared with Clair3's calls (and filling positions Clair3 left empty) in ABO_EVIDENCE.
    -a keeps zero-depth positions in the output so "no coverage" is explicit.
*/
process ABO_PILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)
    path bed

    output:
    tuple val(meta), path("*.pileup.txt"), emit: pileup
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools mpileup \\
        ${args} \\
        -l ${bed} \\
        -f ${fasta} \\
        -o ${prefix}.pileup.txt \\
        ${bam}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileup.txt
    """
}
