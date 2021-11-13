process GATK_SELECT_VARIANTS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    val(variantType)
    val(prefix)
    tuple val(joint_name), path(vcfIndex), path(vcf)
    val(resourceFilesArg)
    path(resourceFiles)
    path(resourceFileIndexes)
    path(reference)
    path("*")

    output:
    tuple val(joint_name), path("*.${prefix}_${variantType}.vcf.gz.tbi"), path("*.${prefix}_${variantType}.vcf.gz"), emit: variantsVcfTuple


    script:

    def excludeIntervalsArg = { resourceFilesArg != "" ? "-XL ${resourceFilesArg} " : "" }

    """
    ${params.gatk_path} SelectVariants --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include ${variantType} \\
        ${excludeIntervalsArg} \\
        ${params.arguments} \\
        -O ${joint_name}.${prefix}_${variantType}.vcf.gz
    """

    stub:

    """
    echo "${params.gatk_path} SelectVariants -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -V ${vcf} \\
        --select-type-to-include ${variantType} \\
        ${finalResourceFilesArg} \\
        ${params.arguments} \\
        -O ${joint_name}.${prefix}_${variantType}.vcf.gz
    """
}
