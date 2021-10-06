
	check_exit $SAMTOOLS stats -F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam -r $REFERENCE > $OUT_DIR/stats/$SAMPLE_ID.recal_reads.SamtoolStats



nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbbwa"
params.save_mode = 'copy'
params.should_publish = true



process process_name {
    tag "something"
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(somefile)

    output:
    path("pattern"),  emit: "ch_output"

    script:

    """
    echo "nothing"
    """

    stub:

    """
    echo "nothing on stub"
    """

}