process MergeCountMatrix {
    cpus 1

    tag "Merge_count_matrix"
    publishDir "$params.outdir/merge", pattern: "*.txt", mode: 'copy'

    input:
    path(COUNT)
    path(INPUT)

    output:
    tuple path("${params.project}.merge.count.txt"),  path("${params.project}.merge.tpm.txt"),  path("${params.project}.merge.fpkm.txt")

    script:
    def project = params.project
    def gene_lens = params.references[params.ref].gene_lens
    """
    merge_count_featureCounts.pl ${INPUT} ./ > ${project}.merge.count.txt
    count2TPM_FPKM_v2.pl ${project}.merge.count.txt ${gene_lens} ${project}.merge
    """
}
