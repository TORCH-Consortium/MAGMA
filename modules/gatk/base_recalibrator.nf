nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/base_recalibrator"
params.save_mode = 'copy'
params.should_publish = true

process GATK_BASE_RECALIBRATOR {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    gatk BaseRecalibrator -Xmx${task.memory.giga}G \\
    --known-sites ${dbsnp} \\
    -R ${ref_fasta} \\
	-I $OUT_DIR/mapped/$SAMPLE_ID.dedup_reads.bam \\
    -O $OUT_DIR/mapped/$SAMPLE_ID.recal_data.table
    """

    stub:

    """
    """
}
