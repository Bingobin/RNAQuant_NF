process INDEX_BAM {
    tag "index bam in $ID"

    input:
    tuple val(ID), path(BAM)

    output:
    tuple val(ID), path(BAM), path("${BAM}.bai")

    script:
    """
    samtools index -@ $task.cpus ${BAM}
    """
}
