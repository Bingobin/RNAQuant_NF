process BAM_TO_BIGWIG {
    tag "bam to bigwig in $ID"
    cpus 4
    publishDir "$params.outdir/bigwig", mode: 'copy'

    input:
    tuple val(ID), path(BAM), path(BAI)

    output:
    tuple val(ID), path("${ID}.CPM.bw"), emit: bigwig

    script:
    """
    bamCoverage \
        --bam ${BAM} \
        --outFileName ${ID}.CPM.bw \
        --normalizeUsing CPM \
        --binSize ${params.bigwig_bin_size} \
        --numberOfProcessors $task.cpus
    """
}
