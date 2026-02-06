include { MergeCountMatrix } from '../modules/merge_counts'

workflow MERGE {
    take:
    counts
    input_file

    main:
    merged = MergeCountMatrix(counts.collect(), input_file)

    emit:
    merged
}
