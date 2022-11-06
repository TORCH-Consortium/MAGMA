process SAMPLESHEET_VALIDATION {

    input:
        path(samplesheet)

    output:
        path(samplesheet)

    script:

        """
        samplesheet_validation.py ${samplesheet}
        """
}
