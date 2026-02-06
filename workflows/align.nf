include { ALIGN_HISAT2 } from '../modules/align_hisat2'
include { ALIGN_STAR } from '../modules/align_star'

workflow ALIGN {
    take:
    fastq

    main:
    if (params.aligner == 'hisat2') {
        ch_align = ALIGN_HISAT2(fastq)
    } else {
        ch_align = ALIGN_STAR(fastq)
    }

    emit:
    bam = ch_align.bam
    stat = ch_align.stat
}
