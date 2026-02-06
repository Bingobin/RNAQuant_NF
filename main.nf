#!/usr/bin/env nextflow

params.input = "$projectDir/bin/samplesheet_test.tsv"
params.project = "RNAseq"
params.outdir = "results"
params.aligner = params.aligner ?: 'star'
params.bam_list = params.bam_list ?: null
params.skip_align = params.skip_align ?: false

// Reference selection (hsa/mmu) is defined in `nextflow.config` under `params.references`.
// The script selects the proper reference entry at runtime.

// Resolve selected reference configuration
def refKey = params.get('ref','hsa')
def ref = params.references ? params.references[refKey] : null
if( !ref ) {
    log.error "Reference with key '${refKey}' not found in params.references.\nPlease set params.ref to one of: ${params.references?.keySet() ?: '[]'}"
    System.exit(1)
}
params.ref_data = ref

log.info """
    R N A S E Q - N F   P I P E L I N E
    ===================================
    selected_ref  :  ${refKey}
    gtf           :  ${params.ref_data.gtf}
    hisat_index   :  ${params.ref_data.hisat_index}
    star_index    :  ${params.ref_data.star_index}
    gene_lens     :  ${params.ref_data.gene_lens}
    input         :  ${params.input}
    outdir        :  ${params.outdir}
    aligner       :  ${params.aligner}
    """.stripIndent(true)

include { QC } from './workflows/qc'
include { ALIGN } from './workflows/align'
include { COUNT } from './workflows/count'
include { MERGE } from './workflows/merge'
include { REPORT } from './workflows/report'
include { PREP_BAM } from './workflows/prep_bam'

workflow {
    if (params.skip_align) {
        if (!params.bam_list) {
            log.error "params.skip_align=true requires --bam_list <tsv> with columns: ID, BAM"
            System.exit(1)
        }

        Channel.fromPath(params.bam_list)
        .splitCsv(header:true, sep: '\t')
        .map{["${it.ID}" , file("${it.BAM}")]}
        .set{ch_bam_in}

        ch_bam = PREP_BAM(ch_bam_in).bam
        ch_count = COUNT(ch_bam)
    } else {
        Channel.fromPath(params.input)
        .splitCsv(header:true, sep: '\t')
        .map{["${it.ID}" ,["${it.R1}", "${it.R2}"]]}
        .set{ch_sample}

        ch_qc = QC(ch_sample)
        ch_align = ALIGN(ch_qc.fastq)
        ch_count = COUNT(ch_align.bam)
    }

    def input_file = file(params.input)
    MERGE(ch_count.counts, input_file)

    REPORT(ch_count.summary)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong" )
}
