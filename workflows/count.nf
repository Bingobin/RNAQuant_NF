include { INFER_STRAND } from '../modules/infer_strand'
include { FeatureCounts } from '../modules/featurecounts'

workflow COUNT {
    take:
    bam

    main:
    ch_strand = INFER_STRAND(bam)
    ch_count = FeatureCounts(bam.join(ch_strand.strand))

    emit:
    counts = ch_count.counts
    summary = ch_count.summary
    out = ch_count.out
    strand = ch_strand
}
