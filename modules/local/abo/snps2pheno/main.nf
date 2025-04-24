process ABOSNPS2PHENO {
    tag "COMPILING ABO RESULTS"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'biocontainers/python:3.9--1' }"

    publishDir "${params.outdir}", mode: 'copy'
    
    /*
    changes to python script not processed properly on re-run
    Disable caching for the process to repeat every time
    */
    // cache false
    
    input:
    path(samples_dir)
    path(per_sample_processing)

    output:
    path "ABO_result.txt"        , emit: txt
    path "ABO_result.xlsx"       , emit: xls
    path "ABO_results.log"       , emit: log
    path "final_export.csv"      , emit: csv
    path "versions.yml"          , emit: versions

    script:
    """
    python3 $projectDir/bin/aggregate_abo_reports.py \\
        per_sample_processing 2>&1 | tee ABO_results.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        xlsxwriter: \$(python -c "import xlsxwriter; print(xlsxwriter.__version__)")
        openpyxl: \$(python -c "import openpyxl; print(openpyxl.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch final_export.csv
    touch ABO_results.log
    touch ABO_result.xlsx
    touch ABO_result.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        xlsxwriter: \$(python -c "import xlsxwriter; print(xlsxwriter.__version__)")
        openpyxl: \$(python -c "import openpyxl; print(openpyxl.__version__)")
    END_VERSIONS
    """
}
