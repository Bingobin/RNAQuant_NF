process FASTQC {
    tag "fastqc in $ID"
    publishDir "$params.outdir/fastqc", pattern: "*.{html,zip}", mode: 'copy'

    input:
    tuple val(ID), path(fastq)

    output:
    tuple val(ID), path("*.{html,zip}")

    script:
    """
    mv ${fastq[0]} ${ID}_R1.fastq.gz
    mv ${fastq[1]} ${ID}_R2.fastq.gz
    fastqc -t $task.cpus ${ID}_R1.fastq.gz ${ID}_R2.fastq.gz
    """
}
