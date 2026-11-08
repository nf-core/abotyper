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
        ? 'oras://community.wave.seqera.io/library/wget:1.25.0--6d2130bfb28314a4'
        : 'community.wave.seqera.io/library/wget:1.25.0--452eb4bbcccd1c30'}"

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
