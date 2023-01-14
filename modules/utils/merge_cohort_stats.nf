process UTILS_MERGE_COHORT_STATS {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        approved_samples_tsv
        rejected_samples_tsv
        call_wf_cohort_stats_tsv

    output:
        path("*merged_cohort_stats.tsv")


    shell:
        """
        generate_merged_cohort_stats.py \\
            --relabundance_approved_tsv ${approved_samples_tsv} \\
            --relabundance_rejected_tsv ${rejected_samples_tsv}\\
            --call_wf_cohort_stats_tsv ${call_wf_cohort_stats_tsv}\\
            --output_file ${params.vcf_name}.merged_cohort_stats.tsv
        """
}
