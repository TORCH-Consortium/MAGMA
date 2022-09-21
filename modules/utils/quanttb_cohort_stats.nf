process UTILS_QUANTTB_COHORT_STATS {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")

    output:
        path("approved_samples.quanttb_cohort_stats.tsv"), emit: approved_samples
        path("rejected_samples.quanttb_cohort_stats.tsv"), emit: rejected_samples

    script:
        """
        quanttb_filter.py *csv
        """
}
