#!/usr/bin/env nextflow

params.input = "$projectDir/bin/samplesheet_test.tsv"
params.project = "RNAseq"
params.outdir = "results"

// Detect CSV/TSV separator from first non-empty line.
def detectSep(String path) {
    def f = file(path)
    def line = f.readLines().find { it?.trim() }
    if (!line) return '\t'
    def hasTab = line.contains('\t')
    def hasComma = line.contains(',')
    if (hasTab && !hasComma) return '\t'
    if (hasComma && !hasTab) return ','
    return '\t'
}

// Reference selection (hsa/mmu) is defined in `nextflow.config` under `params.references`.
// The script selects the proper reference entry at runtime.

include { QC } from './workflows/qc'
include { ALIGN } from './workflows/align'
include { COUNT } from './workflows/count'
include { MERGE } from './workflows/merge'
include { REPORT } from './workflows/report'
include { PREP_BAM } from './workflows/prep_bam'
include { BAM_TO_BIGWIG } from './modules/bam_to_bigwig'

workflow {
    main:
    def refKey = params.get('ref','hsa')
    def ref = params.references ? params.references[refKey] : null
    if( !ref ) {
        log.error "Reference with key '${refKey}' not found in params.references.\nPlease set params.ref to one of: ${params.references?.keySet() ?: '[]'}"
        System.exit(1)
    }

    log.info """
        R N A S E Q - N F   P I P E L I N E
        ===================================
        selected_ref  :  ${refKey}
        gtf           :  ${ref.gtf}
        hisat_index   :  ${ref.hisat_index}
        star_index    :  ${ref.star_index}
        gene_lens     :  ${ref.gene_lens}
        input         :  ${params.input}
        outdir        :  ${params.outdir}
        aligner       :  ${params.aligner}
        bigwig        :  ${params.bigwig}
        """.stripIndent(true)

    if (params.skip_align) {
        if (!params.bam_list) {
            log.error "params.skip_align=true requires --bam_list <tsv> with columns: ID, BAM"
            System.exit(1)
        }

        def bam_sep = detectSep(params.bam_list)
        Channel.fromPath(params.bam_list)
        .splitCsv(header:true, sep: bam_sep)
        .map{["${it.ID}" , file("${it.BAM}")]}
        .set{ch_bam_in}

        ch_bam = PREP_BAM(ch_bam_in).bam
        ch_count = COUNT(ch_bam)
    } else {
        def input_sep = detectSep(params.input)
        Channel.fromPath(params.input)
        .splitCsv(header:true, sep: input_sep)
        .map{["${it.ID}" ,["${it.R1}", "${it.R2}"]]}
        .set{ch_sample}

        ch_qc = QC(ch_sample)
        ch_align = ALIGN(ch_qc.fastq)
        ch_bam = ch_align.bam
        ch_count = COUNT(ch_bam)
    }

    if (params.bigwig == true) {
        BAM_TO_BIGWIG(ch_bam)
    }

    def input_file = file(params.input)
    MERGE(ch_count.counts, input_file)

    REPORT(ch_count.summary)

    onComplete:
    log.info ( workflow.success ? "\nDone! See results --> $params.outdir\n" : "Oops.. someting went wrong" )
}
