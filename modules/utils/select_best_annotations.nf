process UTILS_SELECT_BEST_ANNOTATIONS {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        val(joint_name)
        path("annotations_and_tranches_json_files/*")
        path("recal_vcf_files/*")
        path("tranches_files/*")


    output:
        tuple val(joint_name), path("best_annotation_files/*.tbi"), path("best_annotation_files/*.recal.vcf.gz"), emit: bestRecalVcfTuple
        tuple val(joint_name), path("best_annotation_files/*.tranches"), emit: bestTranchesFile


    script:
        """
        mkdir best_annotation_files

        select_best_annotations.py \\
            --input_directory annotations_and_tranches_json_files \\
            --output_directory best_annotation_files
        """
}
