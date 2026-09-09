/*
  MODULE: HAPLOSCAN
  Description:
    Read-level ABO haplotype scanner using pysam.

    Replaces the two-step  SAMTOOLS_MPILEUP → MPILEUPSTATS  with a single
    Python/pysam process that works directly on the aligned BAM.

    Because ONT reads can span an entire amplicon, variants on the same read
    are confirmed on the same physical molecule — enabling phasing without
    any wet-lab changes.

  Outputs:
    *.AlignmentStatistics.tsv  — identical column layout to stats_from_pileup.py;
                                  all downstream scripts are unchanged.
    *.Haplotypes.tsv           — per-read haplotype table.
    ABOReadPolymorphisms.txt   — same polymorphic-position summary as before.
*/

process HAPLOSCAN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pip_pysam:8f917a7d053340c0'
        : 'community.wave.seqera.io/library/pip_pysam:d1234b67d00cd6a5'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta1), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.AlignmentStatistics.tsv"), emit: tsv
    tuple val(meta), path("*.Haplotypes.tsv"),          emit: haplotypes
    path ("ABOReadPolymorphisms.txt"),                   emit: txt
    tuple val("${task.process}"), val('python'),
          eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pysam'),
          eval('python3 -c "import pysam; print(pysam.__version__)"'), emit: versions_pysam, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pysam_haploscan.py \\
        -b ${bam} \\
        -f ${fasta} \\
        -o ${prefix} \\
        -s ABOReadPolymorphisms.txt \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.AlignmentStatistics.tsv
    touch ${prefix}.Haplotypes.tsv
    touch ABOReadPolymorphisms.txt
    """
}
