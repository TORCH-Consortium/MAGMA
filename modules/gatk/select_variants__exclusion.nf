process GATK_SELECT_VARIANTS__EXCLUSION {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(analysisMode)
        tuple val(joint_name), path(filteredVcfIndex), path(filteredVcf)
        path(intervalFile)
        path(reference)
        path("*")


    output:
        tuple val(joint_name), path("*.filtered_${analysisMode}_exc-rRNA.vcf.gz.tbi"), path("*.filtered_${analysisMode}_exc-rRNA.vcf.gz")

    script:

        """
        ${params.gatk_path} SelectVariants --java-options "-Xmx${task.memory.giga}G" \\
            -V ${filteredVcf} \\
            --exclude-intervals ${intervalFile} \\
            ${params.arguments} \\
            -O ${joint_name}.filtered_${analysisMode}_exc-rRNA.vcf.gz

        """

    stub:

        """
        touch ${joint_name}.filtered_SNP_exc-rRNA.vcf.gz
        touch ${joint_name}.filtered_SNP_exc-rRNA.vcf.gz.tbi
        """
}
