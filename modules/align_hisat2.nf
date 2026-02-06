process ALIGN_HISAT2 {
    tag "Hista2 align in $ID"
    publishDir "$params.outdir/align", pattern: "*.xls", mode: 'copy'

    input:
    tuple val(ID), path(fastq)

    output:
    tuple val(ID), path("${ID}.sorted.bam"), path("${ID}.sorted.bam.bai"), emit: bam
    tuple val(ID), path("${ID}.Hisat2Genome.MapReadsStat.xls"), emit: stat

    script:
    """
    hisat2 -p $task.cpus --phred33 --sensitive -I 1 -X 1000 -x ${params.ref_data.hisat_index} -1 ${fastq[0]} -2 ${fastq[1]} 2> ${ID}.Hisat2Genome.MapReadsStat.xls | samtools view -@ $task.cpus -b -S -o ${ID}.bam -
    samtools sort -@ $task.cpus -o ${ID}.sorted.bam ${ID}.bam
    samtools index ${ID}.sorted.bam
    rm ${ID}.bam
    """
}
