process GATK_VARIANT_RECALIBRATOR {
    tag "annotation: ${annotations}"
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
        path("*.pdf")
        path("*${analysisMode}*.command.log"), emit: annotationsLog

    script:

        def finalResourceFilesArg =    (resourceFilesArg  ? "--resource:${resourceFilesArg}" : "")

        def optionalAnnotationPrefix = ""

        if (task.process.split("__").length == 1) {
            optionalAnnotationPrefix = ""
        } else {
            optionalAnnotationPrefix = ".${task.process.split("__")[-1]}"
        }

        """
        ${params.gatk_path} VariantRecalibrator --java-options "-Xmx${task.memory.giga}G" \\
            -R ${reference} \\
            -V ${variantsVcf} \\
            ${finalResourceFilesArg} \\
            ${annotations} \\
            ${params.arguments} \\
            -mode ${analysisMode} \\
            --tranches-file ${joint_name}.${analysisMode}${optionalAnnotationPrefix}.tranches \\
            --rscript-file ${joint_name}.${analysisMode}${optionalAnnotationPrefix}.R \\
            --output ${joint_name}.${analysisMode}${optionalAnnotationPrefix}.recal.vcf.gz \\
            --output-model ${joint_name}.${analysisMode}${optionalAnnotationPrefix}.model \\
            2>${joint_name}.${analysisMode}${optionalAnnotationPrefix}.command.log

        cp ${joint_name}.${analysisMode}${optionalAnnotationPrefix}.command.log .command.log

        """

    stub:

        """
        touch ${joint_name}.${analysisMode}.tranches
        touch ${joint_name}.${analysisMode}.R
        touch ${joint_name}.${analysisMode}.recal.vcf.gz
        touch ${joint_name}.${analysisMode}.mod
        """
}
