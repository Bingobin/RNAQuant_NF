include { MultiQC } from '../modules/multiqc'

workflow REPORT {
    take:
    done

    main:
    report = MultiQC(done)

    emit:
    report
}
