process UTILS_QUANTTB_COHORT_STATS {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")

    output:
        path("*tsv")

    shell:
        '''
        echo -e "SAMPLE\tREFNAME\tTOTSCORE\tRELABUNDANCE\tRELABUNDANCE_THRESHOLD_MET\tDEPTH\tDERIVED_NAME" > !{params.vcf_name}.quanttb_cohort_stats.tsv

        cat *csv >> !{params.vcf_name}.quanttb_cohort_stats.tsv
        '''
}
