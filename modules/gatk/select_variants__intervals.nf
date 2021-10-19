
process GATK_SELECT_VARIANTS__INCLUSION {
    tag "${analysisMode}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(structuralVariantsIndex), path(structuralVariantsVcf)
    path(intervalsFile)

    output:
    path("*.filtered_${analysisMode}_exc-rRNA.vcf.gz")

    script:

    """
    ${params.gatk_path} SelectVariants -Xmx${task.memory.giga}G \\
        -V ${structuralVariantsVcf} \\
        -L ${intervalsFile} \\
        -O ${sampleName}.potentialSV.DRgenes.vcf.gz
    """

    stub:

    """
    touch ${joint_name}.filtered_SNP_exc-rRNA.vcf.gz
    """
}
