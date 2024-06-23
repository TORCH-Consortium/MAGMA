process UTILS_FASTQ_COHORT_VALIDATION {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("fastq_reports/*")
        path(magma_validated_samplesheet)

//    output:
    //        path("*.fastqs.passed.tsv"), emit: passed_fastqs
    //        path("*.fastqs.failed.tsv"), optional: true

    script:

        """
        csvtk concat fastq_reports/* |  csvtk csv2json > merged_fastq_reports.json

        csvtk csv2json ${magma_validated_samplesheet} -k MagmaSampleName > samplesheet.json

        fastq_cohort_validation.py ${params.vcf_name}
        """


    stub:

        """
        touch ${params.vcf_name}.passed.tsv ${params.vcf_name}.failed.tsv
        """

}
