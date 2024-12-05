process AGGREGATE_REPORTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path('*')

    output:
    path "ABO_result.txt", emit: txt
    path "ABO_result.xlsx", emit: xlsx
    path "versions.yml", emit: versions

    script:
    """
    python $projectDir/bin/Aggregate_ABO_reports.py .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
