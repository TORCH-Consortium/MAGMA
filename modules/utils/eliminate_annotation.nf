process UTILS_ELIMINATE_ANNOTATION {
    tag "${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(analysisType)
        path(annotationsFile)

    output:
        path("*annotations.txt"), emit: reducedAnnotationsFile

    script:
        """
        reduce_annotations.py \\
            -i ${annotationsFile} \\
            -o ${params.vcf_name}.${analysisType}.${task.process.tokenize('__')[-1]}.annotations.txt
        """
}
