
        check_exit $SNP_SITES -o $OUT_DIR/fasta/$JOINT_NAME/$JOINT_NAME.95X.variable.IncComplex.fa $OUT_DIR/fasta/$JOINT_NAME/$JOINT_NAME.95X.IncComplex.fa


nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/snpsites"
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