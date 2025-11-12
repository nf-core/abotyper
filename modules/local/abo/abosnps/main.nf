process GETABOSNPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/json5_pandas_python:ce15c1d2f7290f72' :
        'community.wave.seqera.io/library/json5_pandas_python:e3184b0698afebbd' }"

    input:
    tuple val(meta), path(variants_freq), path(coverage), val(exon_n)

    output:
    tuple val(meta), path("*.ABOPhenotype.txt"), emit: phenotype
    tuple val(meta), path("*.log.txt")         , emit: log
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    predict_abo_phenotype.py \\
        -i ${variants_freq} \\
        -o ${prefix}.ABOPhenotype.txt \\
        -c ${coverage} \\
        -e ${exon_n} \\
        2>&1 | tee ${prefix}.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        json5: \$(python3 -c "import json5; print(json5.__version__)")
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
        json5: \$(python3 -c "import json5; print(json5.__version__)")
    END_VERSIONS
    """
}
