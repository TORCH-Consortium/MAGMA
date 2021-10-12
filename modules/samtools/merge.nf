

$SAMTOOLS merge -f $OUT_DIR/mapped/$SAMPLE_ID.sorted_reads.bam $OUT_DIR/mapped_singles/$SAMPLE_ID.* -@ $SAMTOOLS_THREADS

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
