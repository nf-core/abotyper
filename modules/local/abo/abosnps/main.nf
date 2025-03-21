process GETABOSNPS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(variants_freq)

    output:
    tuple val(meta), path("*.ABOPhenotype.txt"), emit: phenotype
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 $projectDir/bin/predict_abo_phenotype.py \\
        -i ${variants_freq} \\
        -o ${prefix}.ABOPhenotype.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch -o ${prefix}.ABOPhenotype.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
