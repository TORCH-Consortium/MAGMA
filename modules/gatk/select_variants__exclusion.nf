process GATK_SELECT_VARIANTS__EXCLUSION {
    tag "${analysisMode}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    val(analysisMode)
    path(filteredVcf)
    path(intervalFile)

    output:
    path("*.filtered_${analysisMode}_exc-rRNA.vcf.gz")

    script:

    """
    ${params.gatk_path} SelectVariants -Xmx${task.memory.giga}G \\
        -V ${filteredVcf} \\
        --exclude-intervals ${intervalFile} \\
        -O ${joint_name}.filtered_${analysisMode}_exc-rRNA.vcf.gz

    """

    stub:

    """
    touch ${joint_name}.filtered_SNP_exc-rRNA.vcf.gz
    """
}
