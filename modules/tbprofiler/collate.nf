
tb-profiler collate --db $RESISTANCE_DB -d $OUT_DIR/resistance/$JOINT_NAME/XBS/results/ -p $OUT_DIR/resistance/$JOINT_NAME/XBS/$JOINT_NAME.XBS.resistance



nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/tbprofiler/collate"
params.save_mode = 'copy'
params.should_publish = true


process TBPROFILER_COLLATE {
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path("results/*")

    output:
    path("tbprofiler*")

    script:
    """
    tb-profiler update_tbdb
    tb-profiler collate
    cp tbprofiler.txt tbprofiler_cohort_report.tsv
    """

    stub:
    """
    touch tbprofiler_cohort_report.tsv
    """
}


