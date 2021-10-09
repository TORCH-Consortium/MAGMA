nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/combine_gvcfs"
params.save_mode = 'copy'
params.should_publish = true

process GATK_COMBINE_GVCFS {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:


    script:

    """
    gatk CombineGVCFs -Xmx${task.memory.giga}G \\
        -R $REFERENCE \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation $GVCFs \\
        -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.combined.vcf.gz
    """

    stub:

    """
    """
}

