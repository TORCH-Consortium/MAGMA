process MULTIQC {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    label 'cpu_4_memory_16'

    input:
    	path(multiqc_config)
        path("*")

    output:
        tuple path("multiqc_data"), path("multiqc_report.html")


    script:
        def config = multiqc_config ? "--config $multiqc_config" : ''

        def prepare_script_options =  !params.skip_merge_analysis ? '--skip_merge_analysis' : ''

        """
        preprocess_multiqc_input.py   \\
        ${prepare_script_options} \\
        --merged_cohort_stats joint.merged_cohort_stats.tsv \\
	--distance_matrix joint.ExDR.ExComplex.snp_dists.tsv

        ${params.multiqc_path} $config .
        """

    stub:
        """
        mkdir multiqc_data

        touch multiqc_report.html
        """
}
