/*
  MODULE: ABO_EVIDENCE
  Description:
    Combines the two independent lines of evidence for the ABO panel positions:
    Clair3 (ABO_CLAIR2METRICS, from the gVCF) and the BAM pileup (ABO_PILEUP).
    Runs for every sample.

      * writes the table PREDICTABOPHENOTYPE scores: Clair3 is the primary line, and
        the BAM pileup fills any position Clair3 had no evidence for (or every position
        if Clair3 produced nothing)
      * compares the two lines at every panel position and writes a per-position
        concordance table
      * writes a one-line status file (evidence source, covered-by-both counts,
        discordant positions overall and at the primary typing markers, primary
        positions filled from the BAM) that aggregate_abo_reports.py turns into
        Typing-sheet columns and REVIEW flags

    See bin/pileup2metrics.py.
*/
process ABO_EVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pyyaml:10a1a757f29eaf02'
        : 'community.wave.seqera.io/library/python_pyyaml:a1e09ed0f4856f89'}"

    input:
    // staged in clair3/ so its name (<id>.AlignmentStatistics.tsv) does not collide with the
    // merged output of the same name: Nextflow never matches input files as outputs
    tuple val(meta), path(pileup), path(clair_metrics, stageAs: "clair3/*")
    path panel_file

    output:
    tuple val(meta), path("*.AlignmentStatistics.tsv"), emit: tsv
    tuple val(meta), path("*.evidence_status.tsv"),     emit: status
    tuple val(meta), path("*.concordance.tsv"),         emit: concordance
    tuple val(meta), path("*.pileup_metrics.tsv"),      emit: pileup_tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pyyaml'), eval('python3 -c "import yaml; print(yaml.__version__)"'), emit: versions_pyyaml, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def clair  = clair_metrics ? "--clair-metrics ${clair_metrics}" : ''
    """
    pileup2metrics.py \\
        --pileup ${pileup} \\
        ${clair} \\
        --panel ${panel_file} \\
        --sample ${meta.id} \\
        --merged-out ${prefix}.AlignmentStatistics.tsv \\
        --pileup-out ${prefix}.pileup_metrics.tsv \\
        --concordance-out ${prefix}.concordance.tsv \\
        --status-out ${prefix}.evidence_status.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.AlignmentStatistics.tsv ${prefix}.evidence_status.tsv \\
          ${prefix}.concordance.tsv ${prefix}.pileup_metrics.tsv
    """
}
