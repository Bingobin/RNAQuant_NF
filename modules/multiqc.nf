process MultiQC {
    cpus 1

    tag "MultiQC"
    publishDir "$params.outdir/reports", pattern: "*", mode: 'copy'

    input:
    path summaries

    output:
    path "*"

    script:
    """
    multiqc -n ${params.project}.reports ${workflow.launchDir}/${params.outdir}/*/* -f
    """

}
