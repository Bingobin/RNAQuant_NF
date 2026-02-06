include { FASTQC } from '../modules/fastqc'
include { FastpFilter } from '../modules/fastp'

workflow QC {
    take:
    samples

    main:
    ch_fastqc = FASTQC(samples)
    ch_clean = FastpFilter(samples)

    emit:
    fastq = ch_clean.fastq
    fastqc = ch_fastqc
    fastp = ch_clean.stat
}
