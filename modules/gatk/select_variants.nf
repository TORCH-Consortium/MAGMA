/*
FIXME: Documentation comments

*/


process GATK_SELECT_VARIANTS {
    tag "type: ${variantType}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    val(variantType)
    val(resourceFilesArg)
    path("*")

    output:

    script:

    def finalResourceFilesArg =    (resourceFilesArg  ? "-XL:${resourceFilesArg}" : "")

    """
    gatk SelectVariants -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -V ${annotatedVcf} \\
        --select-type-to-include ${variantType} \\
        ${finalResourceFilesArg} \\
        ${params.arguments} \\
        -O ${joint_name}.raw_${variantType}.vcf.gz
    """

    stub:

    """
    """
}
