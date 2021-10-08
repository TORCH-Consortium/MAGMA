nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/select_variants"
params.save_mode = 'copy'
params.should_publish = true


// NOTE: Process-2 for excluding intervals

process GATK_SELECT_VARIANTS_INTERVALS {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:

    script:

    """
    gatk SelectVariants -Xmx${task.memory.giga}G \\
    -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_inc-rRNA.vcf.gz \\
    --exclude-intervals $RRNA \\
    -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_exc-rRNA.vcf.gz

    """

    stub:

    """
    """
}
