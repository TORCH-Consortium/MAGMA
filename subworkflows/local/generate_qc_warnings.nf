process UTILS_GENERATE_QC_WARNINGS {

    tag "joint_name: ${params.vcf_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path merged_cohort_stats_tsv
        val lofreq_cutoff

    output:
        path "warnings.txt", emit: warnings_ch

    script:
    """
    generate_cohort_warnings.py \\
        --merged-cohort-stats ${merged_cohort_stats_tsv} \\
        --lofreq-cutoff ${lofreq_cutoff} \\
        --output warnings.txt
    """
}
