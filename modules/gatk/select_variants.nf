process GATK_SELECT_VARIANTS {
    tag "type: ${variantType}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    val(variantType)
    tuple val(joint_name), path(annotatedVcfIndex), path(annotatedVcf)
    path(reference)
    path("*")

    output:
    tuple val(joint_name), path("*.raw_${variantType}.vcf.gz.tbi"), path("*.raw_${variantType}.vcf.gz")


    script:

    """
    ${params.gatk_path} SelectVariants --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -V ${annotatedVcf} \\
        --select-type-to-include ${variantType} \\
        ${params.arguments} \\
        -O ${joint_name}.raw_${variantType}.vcf.gz
    """

    stub:

    """
    echo "${params.gatk_path} SelectVariants -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -V ${annotatedVcf} \\
        --select-type-to-include ${variantType} \\
        ${finalResourceFilesArg} \\
        ${params.arguments} \\
        -O ${joint_name}.raw_${variantType}.vcf.gz"
    """
}
