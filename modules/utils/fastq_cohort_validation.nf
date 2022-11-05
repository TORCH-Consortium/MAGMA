process UTILS_FASTQ_COHORT_VALIDATION {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")

    output:
        path("*.check.passed.tsv"), emit: passed_fastqs
        path("*.check.failed.tsv"), emit: failed_fastqs

    shell:
       
        '''
        cat *passed* > !{params.vcf_name}.fastqs.passed.tsv
        cat *failed* > !{params.vcf_name}.fastqs.failed.tsv
        '''

    stub: 

        """
        touch ${params.vcf_name}.passed.tsv ${params.vcf_name}.failed.tsv
        """ 

}
