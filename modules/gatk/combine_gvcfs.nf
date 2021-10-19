process GATK_COMBINE_GVCFS {
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    tuple val(joint_name), path(gvcfs)
    path(ref_fasta)

    output:
    tuple val(joint_name), path("*.combined.vcf.gz")


    script:

    """
    ${params.gatk_path} CombineGVCFs -Xmx${task.memory.giga}G \\
        -R ${ref_fasta} \\
        ${params.arguments} \\
        ${gvcfs} \\
        -O ${joint_name}.combined.vcf.gz
    """

    stub:

    """
    touch ${joint_name}.combined.vcf.gz
    """
}

