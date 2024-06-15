process SAMPLESHEET_VALIDATION {

    input:
        path(samplesheet)

    output:
        val true

    script:

        """
        samplesheet_validation.py ${samplesheet}
        """
}
