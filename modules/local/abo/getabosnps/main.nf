process ABO_GETABOSNPS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8d69e246c0a530fa88ba496bf8a62bd282e770fb4b68d2877ab777ce4943fed1/data'
        : 'community.wave.seqera.io/library/json5_pandas_python:e3184b0698afebbd'}"

    input:
    tuple val(meta), path(variants_freq), path(coverage)

    output:
    tuple val(meta), path("*.ABOPhenotype.txt"), emit: phenotype
    tuple val(meta), path("*.log.txt"), emit: log
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('pandas'), eval('python3 -c "import pandas; print(pandas.__version__)"'), emit: versions_pandas, topic: versions
    tuple val("${task.process}"), val('json5'), eval('python3 -c "import json5; print(json5.__version__)"'), emit: versions_json5, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    predict_abo_phenotype.py \\
        -i ${variants_freq} \\
        -o ${prefix}.ABOPhenotype.txt \\
        -c ${coverage} \\
        ${args} \\
        2>&1 | tee ${prefix}.log.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ABOPhenotype.txt
    touch ${prefix}.log.txt
    """
}
