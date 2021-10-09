nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/genotype_gvcfs"
params.save_mode = 'copy'
params.should_publish = true

process GATK_GENOTYPE_GVCFS {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:


    output:


    script:

    """
    gatk GenotypeGVCFs -Xmx${task.memory.giga}G \\
        -R $REFERENCE \\
        -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.combined.vcf.gz \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        --sample-ploidy 1 \\
        -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf.gz
    """

    stub:

    """
    """
}

