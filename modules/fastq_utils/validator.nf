process FASTQ_VALIDATOR {
    tag "${sampleName}"
    /* publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish */

    input:
        tuple val(sampleName), path(reads)

    /* output: */
    /*     tuple val(sampleName), path("*.potentialSV.vcf.gz") */

    script:

        """
        ${params.fastq_validator_path} reads[0] reads[1] 
        """

    /* stub: */

    /*     """ */
    /*     touch ${sampleName}.potentialSV.vcf.gz */
    /*     """ */

}
