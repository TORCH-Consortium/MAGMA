nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbbwa"
params.save_mode = 'copy'
params.should_publish = true

process process_name {
    tag "something"
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    $IQTREE -s $OUT_DIR/fasta/$JOINT_NAME/$JOINT_NAME.95X.variable.IncComplex.fa \\
    -T $IQTREE_THREADS \\
    -allnni \\
    --prefix $OUT_DIR/phylogeny/$JOINT_NAME/$JOINT_NAME.95X.IncComplex \\
    #-bb 10000 \\
    -redo
    """

    stub:

    """
    """

}
