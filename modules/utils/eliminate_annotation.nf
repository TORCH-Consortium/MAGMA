process UTILS_ELIMINATE_ANNOTATION {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(joint_name)
        val(analysisType)
        path(annotationsLog)
        path(tranchesFile)

    output:
        path("*annotations.txt"), emit: reducedAnnotationsFile
        tuple val(joint_name), path("*annotations_tranches.json"), emit: annotationsTranchesFile

    script:

        def annotationPrefix = "${task.process.split("__")[-1]}"

        """
        reduce_annotations.py \\
            --input ${annotationsLog} \\
            --output ${joint_name}.${analysisType}.${annotationPrefix}.annotations.txt \\
            --tranches ${tranchesFile} \\
            --output_tranches_and_annotations ${joint_name}.${analysisType}.${annotationPrefix}.annotations_tranches.json
        """
}
