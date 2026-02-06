process MergeCountMatrix {
    tag "Merge_count_matrix"
    publishDir "$params.outdir/merge", pattern: "*.txt", mode: 'copy'

    input:
    path(COUNT)
    path(INPUT)

    output:
    tuple path("${params.project}.merge.count.txt"),  path("${params.project}.merge.tpm.txt"),  path("${params.project}.merge.fpkm.txt")

    script:
    """
    merge_count_featureCounts.pl ${INPUT} ./ > ${params.project}.merge.count.txt
    sed  -i 's/\.[0-9]\+\t/\t/' ${params.project}.merge.count.txt
    grep -v "PAR" ${params.project}.merge.count.txt > a && mv a  ${params.project}.merge.count.txt
    count2TPM_FPKM_v2.pl ${params.project}.merge.count.txt ${params.ref_data.gene_lens} ${params.project}.merge
    """
}
