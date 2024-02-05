process UTILS_FASTQ_COHORT_VALIDATION {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    storeDir params.store_dir

    input:
        path("*")

    output:
        path("*.fastqs.passed.tsv"), emit: passed_fastqs
        path("*.fastqs.failed.tsv"), optional: true

    shell:
       
        '''

        if ls *check.passed* 1> .null 2>&1; then
            cat *check.passed* > !{params.vcf_name}.fastqs.passed.tsv
        else
            echo "No samples passed!"
        fi


        # Creating this file anyhow, since this is an optional output
        touch !{params.vcf_name}.fastqs.failed.tsv

        if ls *check.failed* 1> .null 2>&1; then
            cat *check.failed* > !{params.vcf_name}.fastqs.failed.tsv
        fi
        '''

    stub: 

        """
        touch ${params.vcf_name}.passed.tsv ${params.vcf_name}.failed.tsv
        """ 

}
