process UTILS_ELIMINATE_ANNOTATION {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path(annotationsFile)

    output:
        path("reduced_ordered_annotations.txt")

    script:
        """
        reduce_annotations.py -i ${annotationsFile}
        """
}
