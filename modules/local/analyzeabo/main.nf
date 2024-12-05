process ANALYZE_ABO {
    tag "$meta.id"
    label 'process_medium'

    conda "assets/conda.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR_CONTAINER_HERE' :
        'YOUR_DOCKER_CONTAINER_HERE' }"

    input:
    tuple val(meta), path(reads)
    path reference_exon6
    path reference_exon7
    path alleles

    output:
    tuple val(meta), path("${meta.id}/*/ABOPhenotype.txt"), emit: phenotype
    tuple val(meta), path("${meta.id}/*/ABOReadPolymorphisms.txt"), emit: polymorphisms
    tuple val(meta), path("${meta.id}/*/ReadAlignmentSpreadsheet.csv"), emit: alignment_spreadsheet
    path "${meta.id}/*/alignment/*", emit: alignment_files
    path "${meta.id}_*.log.txt", emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python $projectDir/bin/AnalyzeAbo_Main.py \\
        --reference=$reference_exon6 \\
        --alleles=$alleles \\
        --output=${meta.id}/exon6 \\
        --analysis-type=READS \\
        --reads=$reads

    python $projectDir/bin/AnalyzeAbo_Main.py \\
        --reference=$reference_exon7 \\
        --alleles=$alleles \\
        --output=${meta.id}/exon7 \\
        --analysis-type=READS \\
        --reads=$reads

    mv ${meta.id}/exon6/log.txt ${meta.id}_exon6.log.txt
    mv ${meta.id}/exon7/log.txt ${meta.id}_exon7.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}