process UTILS_ELIMINATE_ANNOTATION {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(analysisType)
        path(annotationsFile)
        path(tranchesFile)

    output:
        path("*annotations.txt"), emit: reducedAnnotationsFile
        path("*annotations_tranches.txt"), emit: annotationsTranchesFile

    script:
        """
        reduce_annotations.py \\
            --input ${annotationsFile} \\
            --output ${params.vcf_name}.${analysisType}.${task.process.tokenize('__')[-1]}.annotations.txt \\
            --tranches ${tranchesFile} \\
            --output_tranches_and_annotations ${params.vcf_name}.${analysisType}.${task.process.tokenize('__')[-1]}.annotations_tranches.txt
        """
}
