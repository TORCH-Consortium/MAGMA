process MULTIQC {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    label 'cpu_low_memory_high'

    input:
        path("*")

    output:
        tuple path("multiqc_data"), path("multiqc_report.html")


    script:

        """
        ${params.multiqc_path} .
        """

    stub:
        """
        mkdir multiqc_data

        touch multiqc_report.html
        """
}
