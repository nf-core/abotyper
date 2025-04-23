process MAKEINDEX {
    tag "FAI to BED"
    label 'process_single'

    input:
    tuple val(meta), path(exon6fai)
    tuple val(meta1), path(exon7fai)

    output:
    tuple val(meta), path("*_exon6.bed"),  emit: exon6bed
    tuple val(meta1), path("*_exon7.bed"), emit: exon7bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id.toString().tokenize('.')[0]}"
    def prefix1 = task.ext.prefix1 ?: "${meta1.id.toString().tokenize('.')[0]}"

    """
    awk -v FS="\t" -v OFS="\t" '{print \$1 FS "0" FS (\$2)-1}' $exon6fai > ${prefix}_exon6.bed
    awk -v FS="\t" -v OFS="\t" '{print \$1 FS "0" FS (\$2)-1}' $exon7fai > ${prefix1}_exon7.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id.toString().tokenize('.')[0]}"
    def prefix1 = task.ext.prefix1 ?: "${meta1.id.toString().tokenize('.')[0]}"

    """
    touch ${prefix}_exon6.bed
    touch ${prefix1}_exon7.bed
    """
}