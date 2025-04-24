process GETABOSNPS {
    tag "$meta.id"
    label 'process_single'
    
    /*
    changes to python script not processed properly on re-run
    Disable caching for the process to repeat every time
    */
    // cache false  

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.3.5' :
        'biocontainers/pandas:1.3.5' }"

    input:
    tuple val(meta), path(variants_freq), path(coverage)
    val exon_n

    output:
    tuple val(meta), path("*.ABOPhenotype.txt"), emit: phenotype
    tuple val(meta), path("*.log.txt")         , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 $projectDir/bin/predict_abo_phenotype.py \\
        -i ${variants_freq} \\
        -o ${prefix}.ABOPhenotype.txt \\
        -c ${coverage} \\
        -e ${exon_n} \\
        2>&1 | tee ${prefix}.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ABOPhenotype.txt
    touch ${prefix}.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
