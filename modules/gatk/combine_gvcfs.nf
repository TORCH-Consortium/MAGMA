process GATK_COMBINE_GVCFS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    val(joint_name)
    val(gvcfs_string)
    path(gvcfs)
    path(ref_fasta)
    path("*")

    output:
    tuple val(joint_name),  path("*.combined.vcf.gz.tbi"), path("*.combined.vcf.gz")


    script:

    """
    ${params.gatk_path} CombineGVCFs --java-options "-Xmx${task.memory.giga}G" \\
        -R ${ref_fasta} \\
        ${params.arguments} \\
        --variant ${gvcfs_string} \\
        -O ${joint_name}.combined.vcf.gz
    """

    stub:

    """
    touch ${joint_name}.combined.vcf.gz
    """
}

