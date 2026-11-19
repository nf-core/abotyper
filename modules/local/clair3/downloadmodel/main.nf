/*
  MODULE: CLAIR3_DOWNLOADMODEL
  Description:
    Downloads the two Clair3 PyTorch model files (pileup.pt, full_alignment.pt)
    from a model URL into one directory, passed to CLAIR3 as its user_model.
    Runs once per run; the directory is broadcast to every sample.
*/
process CLAIR3_DOWNLOADMODEL {
    tag { url.tokenize('/').last() }
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a16db78e0a9f742ece8c129ab8ee925e416d02df054fe023c44c32fbaede944c/data'
        : 'community.wave.seqera.io/library/wget:1.25.0--5aec898b4eacb870'}"

    input:
    val url

    output:
    path "clair3_model_dir", emit: model_dir
    tuple val("${task.process}"), val('wget'), eval("wget --version | sed '1!d; s/^GNU Wget //; s/ .*//'"), emit: versions_wget, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def base_url = url.toString().replaceAll('/+$', '')
    """
    mkdir -p clair3_model_dir

    wget ${args} "${base_url}/pileup.pt"         -O clair3_model_dir/pileup.pt
    wget ${args} "${base_url}/full_alignment.pt" -O clair3_model_dir/full_alignment.pt
    """

    stub:
    """
    mkdir -p clair3_model_dir
    touch clair3_model_dir/pileup.pt clair3_model_dir/full_alignment.pt
    """
}
