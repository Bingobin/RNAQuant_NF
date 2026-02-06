process FastpFilter {
    tag "fastp in $ID"
    publishDir "$params.outdir/clean", pattern: "*.{html,json}", mode: 'copy'

    input:
    tuple val(ID), path(fastq)

    output:
    tuple val(ID), path("${ID}_clean_R*.fastq.gz"), emit: fastq
    tuple val(ID), path("*.{html,json}"), emit: stat

    script:
    """
    mv ${fastq[0]} ${ID}_R1.fastq.gz
    mv ${fastq[1]} ${ID}_R2.fastq.gz
    fastp -i ${ID}_R1.fastq.gz -I  ${ID}_R2.fastq.gz -o ${ID}_clean_R1.fastq.gz -O ${ID}_clean_R2.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 -q 20 -u 40 -n 5  -l 50 -w $task.cpus -h ${ID}_fastp.html -j ${ID}_fastp.json
    """

}
