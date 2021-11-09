process GATK_APPLY_VQSR {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    val(analysisMode)
    tuple val(joint_name), path(variantsVcfIndex), path(variantsVcf), path(recalVcfIndex), path(recalVcf)
    path(reference)
    path("*")

    output:
    path("*.filtered_${analysisMode}_inc-rRNA.vcf.gz")

    script:

    """
    ${params.gatk_path} ApplyVQSR --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -V ${variantsVcf} \\
        --tranches-file ${joint_name}.${analysisMode}.tranches \\
        --recal-file ${recalVcf} \\
        ${params.arguments} \\
        -mode ${analysisMode} \\
        -O ${joint_name}.filtered_${analysisMode}_inc-rRNA.vcf.gz
    """

    stub:

    """
    touch ${joint_name}.filtered_${analysisMode}_inc-rRNA.vcf.gz
    """
}

