process ALIGN_STAR {
    tag "STAR  align in $ID"
    publishDir "$params.outdir/align", pattern: "${ID}_{Log.*,SJ.*}", mode: 'copy'
    publishDir "$params.outdir/align", pattern: "${ID}_Aligned.*", mode: 'copy'

    input:
    tuple val(ID), path(fastq)

    output:
    tuple val(ID), path("${ID}_Aligned.sortedByCoord.out.bam"), path("${ID}_Aligned.sortedByCoord.out.bam.bai"), emit: bam
    tuple val(ID), path("${ID}_Log.*"), path("${ID}_SJ.*"), emit: stat

    script:
    """
    STAR --runThreadN $task.cpus --genomeDir ${params.ref_data.star_index} --readFilesIn ${fastq[0]} ${fastq[1]} --readFilesCommand zcat --outFileNamePrefix ${ID}_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
    samtools index ${ID}_Aligned.sortedByCoord.out.bam
    """

}
