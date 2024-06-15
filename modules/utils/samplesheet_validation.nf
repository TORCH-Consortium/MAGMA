process SAMPLESHEET_VALIDATION {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path(samplesheet)

    output:
        val true

    script:

        """
        samplesheet_validation.py ${samplesheet}
        """
}
