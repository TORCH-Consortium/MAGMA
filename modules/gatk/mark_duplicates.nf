nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/mark_duplicates"
params.save_mode = 'copy'
params.should_publish = true


process GATK_MARK_DUPLICATES {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    gatk MarkDuplicates -Xmx${task.memory.giga}G \\
    --METRICS_FILE $OUT_DIR/stats/$SAMPLE_ID.MarkDupMetrics \\
    -I $OUT_DIR/mapped/$SAMPLE_ID.sorted_reads.bam \\
    -O $OUT_DIR/mapped/$SAMPLE_ID.dedup_reads.bam
    """

    stub:

    """
    """

}
