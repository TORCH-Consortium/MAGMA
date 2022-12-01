process TBPROFILER_COLLATE {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        val(joint_name)
        path("results/*")
        path(resistanceDb)

    output:
        path("*${params.prefix}*"), emit: cohort_results
        path("results"), emit: per_sample_results

    script:
        def optionalDb  = resistanceDb ? "--db ${resistanceDb}" : ""

        """

        ${params.tbprofiler_path} collate \\
            ${optionalDb} \\
            -p ${joint_name}.${params.prefix}
        """

    stub:
        """
        touch ${joint_name}.${params.prefix}.txt
        """
}


