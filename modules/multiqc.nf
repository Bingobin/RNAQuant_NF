process MultiQC {
    tag "MultiQC"
    publishDir "$params.outdir/reports", pattern: "*", mode: 'copy'

    input:
    val done

    output:
    path "*"

    script:
    """
    #multiqc -n ${params.project}.fastqc.reports ${workflow.launchDir}/${params.outdir}/fastqc/* -f
    #multiqc -n ${params.project}.clean.reports ${workflow.launchDir}/${params.outdir}/clean/* -f
    #multiqc -n ${params.project}.align.reports ${workflow.launchDir}/${params.outdir}/align/* -f
    #multiqc -n ${params.project}.count.reports ${workflow.launchDir}/${params.outdir}/count/* -f
    multiqc -n ${params.project}.reports ${workflow.launchDir}/${params.outdir}/*/* -f
    """

}
