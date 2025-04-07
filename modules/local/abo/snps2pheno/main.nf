process ABOSNPS2PHENO {
    tag "COMPILING ABO RESULTS"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'biocontainers/python:3.9--1' }"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(sample_dir), val(metas), path(files)

    output:
    path "*.txt", emit: text_files
    path "*.xlsx", emit: excel_files
    path "*.csv", emit: csv_files
    path "*.log", emit: log_files
    path "final_export.csv", emit: final_export
    path "versions.yml", emit: versions

    script:
    def sample_id = sample_dir.name
    """
    mkdir -p ${sample_id}/exon6 ${sample_id}/exon7
    
    for i in {0..1}; do
        if [[ "\${metas[\$i].exon}" == "exon6" ]]; then
            cp ${files[i]} ${sample_id}/exon6/${sample_id}.ABOPhenotype.txt
        elif [[ "\${metas[\$i].exon}" == "exon7" ]]; then
            cp ${files[i]} ${sample_id}/exon7/${sample_id}.ABOPhenotype.txt
        fi
    done

    python3 $projectDir/bin/aggregate_abo_reports.py \\
        ${sample_id} 2>&1 | tee ABO_results.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch final_export.csv
    touch ABO_results.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
