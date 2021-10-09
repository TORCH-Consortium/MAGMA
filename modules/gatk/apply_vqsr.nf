nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/apply_vqsr"
params.save_mode = 'copy'
params.should_publish = true


process GATK_APPLY_VQSR {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

    output:

    script:

    """
    gatk ApplyVQSR -Xmx${task.memory.giga}G \\
        -R $REFERENCE \\
        -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz \\
        --tranches-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.tranches \\
        --recal-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.recal.vcf.gz \\
        --ts-filter-level 99.90 \\
        -AS \\
        --exclude-filtered -mode SNP \\
        -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_inc-rRNA.vcf.gz
    """

    stub:

    """
    """
}

