nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/quanttb"
params.save_mode = 'copy'
params.should_publish = true



process process_name {
    tag ""
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    quanttb quant $FASTQ_FILES -o $OUT_DIR/quanttb/$SAMPLE_ID. -k
    """

    stub:

    """
    """

}
