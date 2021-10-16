/*
FIXME: Documentation comments

*/



process GATK_VARIANT_RECALIBRATOR {
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    val(analysisMode)
    path(variantsVcf)
    val(resourceFilesArg)
    path("*")

    output:
    path("*.recal.vcf.gz"), emit: recalVcf
    path("*.tranches")
    path("*.R")
    path("*.model")

    script:

    def finalResourceFilesArg =    (resourceFilesArg  ? "--resource:${resourceFilesArg}" : "")

    """
    gatk VariantRecalibrator -Xmx${task.memory.giga}G \\
        -R ${reference} \\
        -V ${variantsVcf} \\
        ${finalResourceFilesArg} \\
        ${params.arguments} \\
        -mode ${analysisMode} \\
        --tranches-file ${joint_name}.${analysisMode}.tranches \\
        --rscript-file ${joint_name}.${analysisMode}.R \\
        --output ${joint_name}.${analysisMode}.recal.vcf.gz \\
        --output-model ${joint_name}.${analysisMode}.model
    """

    stub:

    """
    """
}
