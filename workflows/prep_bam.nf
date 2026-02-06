include { INDEX_BAM } from '../modules/index_bam'

workflow PREP_BAM {
    take:
    bam

    main:
    indexed = INDEX_BAM(bam)

    emit:
    bam = indexed
}
