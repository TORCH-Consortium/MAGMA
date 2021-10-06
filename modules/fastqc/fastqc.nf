$FASTQC $SEQ_R1 -t $FASTQC_THREADS -o $OUT_DIR/fastqc/


nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/fastqc"
params.save_mode = 'copy'
params.should_publish = true

process FASTQC {
    tag "${genomeName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(genomeName), path(genomeReads)

    output:
    tuple path('*.html'), path('*.zip')


    script:

    """
    fastqc *fastq*
    """

    stub:
    """
    touch ${genomeName}.html

    touch ${genomeName}.zip
    """
}
