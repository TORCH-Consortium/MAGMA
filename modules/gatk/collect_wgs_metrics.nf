nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/haplotype_caller"
params.save_mode = 'copy'
params.should_publish = true
params.arguments = " --READ_LENGTH 0 --COVERAGE_CAP 10000 --COUNT_UNPAIRED"

process GATK_COLLECT_WGS_METRICS {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    gatk CollectWgsMetrics -Xmx${task.memory.giga}G \\
    -R $REFERENCE \\
    -I $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam \\
    ${params.arguments} \\
    -O $OUT_DIR/stats/$SAMPLE_ID.recal_reads.WgsMetrics
    """

    stub:

    """
    """
}

