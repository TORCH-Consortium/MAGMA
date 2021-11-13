process GATK_MERGE_VCFS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(snpVcfIndex), path(snpVcf), path(indelVcfIndex), path(indelVcf)

    output:
    tuple val(joint_name), path("*.filtered_SNP.RawIndels.vcf.gz.tbi"), path("*.filtered_SNP.RawIndels.vcf.gz")

    script:

    """
    ${params.gatk_path} MergeVcfs --java-options "-Xmx${task.memory.giga}G" \\
        -I ${snpVcf} \\
        -I ${indelVcf} \\
        -O ${joint_name}.filtered_SNP.RawIndels.vcf.gz
    """

    stub:

    """

    """
}

