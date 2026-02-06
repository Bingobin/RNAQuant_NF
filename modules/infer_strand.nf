process INFER_STRAND {
    cpus 2

    tag "infer_experiment.py in $ID"
    publishDir "$params.outdir/count", pattern: "${ID}.strand.result.txt", mode: 'copy'

    input:
    tuple val(ID), path(BAM), path(BAI)

    output:
    tuple val(ID), stdout, emit: strand
    path("${ID}.strand.result.txt"), emit: result

    script:
    """
    infer_experiment.py -i $BAM -r ${params.ref_data.gene_bed} > ${ID}.strand.result.txt
    FWD=\$(cat ${ID}.strand.result.txt | grep "1++,1--,2+-,2-+" | cut -d':' -f2 | tr -d ' ')
    REV=\$(cat ${ID}.strand.result.txt | grep "1+-,1-+,2++,2--" | cut -d':' -f2 | tr -d ' ')
    echo | awk -v f=\$FWD -v r=\$REV '{if (r > 0.6) printf "2"; else if (f > 0.6) printf "1"; else printf "0";}'
    """
}
