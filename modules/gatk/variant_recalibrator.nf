process GATK_VARIANT_RECALIBRATOR {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
    val(analysisMode)
    val(annotations)
    tuple val(joint_name), path(variantsIndex), path(variantsVcf)
    val(resourceFilesArg)
    path(resourceFiles)
    path(resourceFileIndexes)
    path(reference)
    path("*")

    output:
    tuple val(joint_name), path("*.tbi"), path("*.recal.vcf.gz"), emit: recalVcfTuple
    tuple val(joint_name), path("*.tranches"), emit: tranchesFile
    path("*.R")
    path("*.model")

    script:

    def finalResourceFilesArg =    (resourceFilesArg  ? "--resource:${resourceFilesArg}" : "")

    """
    ${params.gatk_path} VariantRecalibrator --java-options "-Xmx${task.memory.giga}G" \\
        -R ${reference} \\
        -V ${variantsVcf} \\
        ${finalResourceFilesArg} \\
        ${annotations} \\
        ${params.arguments} \\
        -mode ${analysisMode} \\
        --tranches-file ${joint_name}.${analysisMode}.tranches \\
        --rscript-file ${joint_name}.${analysisMode}.R \\
        --output ${joint_name}.${analysisMode}.recal.vcf.gz \\
        --output-model ${joint_name}.${analysisMode}.model
    """

    stub:

    """
    touch ${joint_name}.${analysisMode}.tranches
    touch ${joint_name}.${analysisMode}.R
    touch ${joint_name}.${analysisMode}.recal.vcf.gz
    touch ${joint_name}.${analysisMode}.mod
    """
}
