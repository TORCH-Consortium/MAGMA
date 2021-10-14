/*
FIXME: Documentation comments

*/


process GATK_SELECT_VARIANTS {
    tag "${params.type}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:

    output:

    script:


    """
    gatk SelectVariants -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -V ${annotatedVcf} \\
        --select-type-to-include ${params.type} \\
        ${params.arguments} \\
        -O ${joint_name}.raw_${params.type}.vcf.gz
    """

    stub:

    """
    """
}
