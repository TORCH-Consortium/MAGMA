nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/flag_stat"
params.save_mode = 'copy'
params.should_publish = true

process GATK_FLAG_STAT {
    tag ""
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    gatk FlagStat -Xmx${task.memory.giga}G \\
        -R $REFERENCE \\
        -I $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam  \\
        > $OUT_DIR/stats/$SAMPLE_ID.recal_reads.FlagStat

    """

    stub:

    """
    """
}

