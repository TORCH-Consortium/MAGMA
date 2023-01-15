process SAMPLESHEET_VALIDATION {
//TODO: Test with relative paths for input samplesheet

    input:
        path(samplesheet)

    output:
        path(samplesheet)

    script:

        """
        samplesheet_validation.py ${samplesheet}
        """
}
