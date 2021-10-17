
process GATK_SELECT_VARIANTS__INCLUSION {
    tag "${analysisMode}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(potentialVariantsVcf)
    path(intervalsFile)

    output:
    path("*.filtered_${analysisMode}_exc-rRNA.vcf.gz")

    script:

    """


//-V PREFIX.potentialSV.vcf.gz -O PREFIX.potentialSV.DRgenes.vcf.gz -L $RESOURCE_PATH/DRgenes.list

    gatk SelectVariants -Xmx${task.memory.giga}G \\
        -V ${potentialVariantsVcf} \\
        -L ${intervalsFile} \\
        -O ${joint_name}.filtered_${analysisMode}_exc-rRNA.vcf.gz

    """

    stub:

    """
    touch ${joint_name}.filtered_SNP_exc-rRNA.vcf.gz
    """
}
