process NTMPROFILER_COLLATE {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(joint_name)
        path("results/*")

    output:
        path("*${params.prefix}*"), emit: cohort_results

    script:

        """

        ${params.ntmprofiler_path} collate \\
            -d results \\
            -o ${joint_name}.${params.prefix}.txt
        """

    stub:
        """
        touch ${joint_name}.${params.prefix}.txt
        """
}


