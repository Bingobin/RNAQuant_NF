process FeatureCounts {
    tag "FeatrueCounts in $ID"
    publishDir "$params.outdir/count", pattern: "${ID}.*", mode: 'copy'

    input:
    tuple val(ID), path(BAM), path(BAI), val(strand)

    output:
    path("${ID}.count"), emit: counts
    path("${ID}.count.summary"), emit: summary
    path("${ID}.count.out"), emit: out

    script:
    """
    featureCounts -T $task.cpus -a ${params.ref_data.gtf} -o ${ID}.count -p  --countReadPairs --verbose -g gene_id ${BAM} -s ${strand} -B -C  -t exon 2> ${ID}.count.out
    """
}
