process UTILS_FASTQ_COHORT_VALIDATION {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")

    output:
        path("*.fastqs.passed.tsv"), emit: passed_fastqs
        path("*.fastqs.failed.tsv"), optional: true

    script:

        """
        fastq_cohort_validation.py
        """


    stub:

        """
        touch ${params.vcf_name}.passed.tsv ${params.vcf_name}.failed.tsv
        """

}
