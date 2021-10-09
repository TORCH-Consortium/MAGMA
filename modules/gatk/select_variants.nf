nextflow.enable.dsl = 2


params.results_dir = "${params.outdir}/gatk4/select_variants"
params.save_mode = 'copy'
params.should_publish = true


// NOTE: Process-1 for SNP/INDELS

process GATK_SELECT_VARIANTS {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:

    script:


    //TODO: Accomodate in design the type-to-include and file name depends on the variant type SNP/INDEL
    """
    gatk SelectVariants -Xmx${task.memory.giga}G \\
        -R $REFERENCE \\
        -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf.gz \\
        --select-type-to-include SNP \\
        --remove-unused-alternates \\
        --exclude-non-variants \\
        -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz

    gatk SelectVariants -Xmx${task.memory.giga}G \\
        -R $REFERENCE \\
        -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf.gz \\
        --select-type-to-include INDEL \\
        --remove-unused-alternates \\
        --exclude-non-variants \\
        -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_INDEL.vcf.gz

    """

    stub:

    """
    """
}
