process SAMPLESHEET_VALIDATION {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    //NOTE: If this process fails, terminate the pipeline execution immediately
    errorStrategy = 'terminate'

    input:
        path(samplesheet)

    output:
        val true, emit: status
        path("samplesheet.format_valid.csv"), emit: validated_samplesheet

    script:

        """
        samplesheet_validation.py ${samplesheet} samplesheet.format_valid.csv
        """
}
