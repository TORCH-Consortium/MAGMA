process GATK_INDEX_FEATURE_FILE__COHORT {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    tuple val(joint_name), path(genotypedVcf)


    output:
    tuple val(joint_name), path("*.vcf.gz.tbi"), path(genotypedVcf)


    script:

    """
    ${params.gatk_path} IndexFeatureFile --java-options "-Xmx${task.memory.giga}G" \\
        -I ${genotypedVcf}
    """

    stub:

    """
    touch ${joint_name}.vcf.gz
    touch ${joint_name}.vcf.gz.tbi

    """
}
