nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/gatk4/index_feature_file"
params.save_mode = 'copy'
params.should_publish = true


process GATK_INDEX_FEATURE_FILE {
    tag ""

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:


    output:


    script:

    """
    gatk IndexFeatureFile -Xmx${task.memory.giga}G \\
    -I $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf.gz
    """

    stub:

    """

    """
}
