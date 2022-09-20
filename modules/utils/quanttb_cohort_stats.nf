process UTILS_QUANTTB_COHORT_STATS {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("*")

    output:
        path("*tsv")

    script:
        """
        quanttb_filter.py *csv
        """
}
