/*
 * Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
 *
 * This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
 *
 * For quick overview of GPL-3 license, please refer
 * https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
 *
 * - You MUST keep this license with original authors in your copy
 * - You MUST acknowledge the original source of this software
 * - You MUST state significant changes made to the original software
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program . If not, see <http://www.gnu.org/licenses/>.
 */
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
        tuple val(joint_name), path("*${analysisMode}*.command.log"), emit: annotationsLog

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
