process UTILS_FILTER_LOFREQ_TBPROFILER_JSONS_BY_COVERAGE {

    input:
        path resistance_jsons
        path call_wf_cohort_stats_tsv
        val cutoff

    output:
        path "filtered_results/*.json", emit: filtered_jsons

    script:
    """
    python3 ${projectDir}/bin/filter_lofreq_tbprofiler_jsons_by_coverage.py \\
        --resistance-jsons ${resistance_jsons} \\
        --cohort-stats ${cohort_stats_tsv} \\
        --cutoff ${cutoff} \\
        --outdir filtered_results
    """
}
