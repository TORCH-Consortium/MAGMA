process UTILS_SAMPLE_STATS {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(samtoolsStats), path(wgsMetrics), path(flagStats), path(ntmFraction)

    output:
        path("*.stats.tsv")

    script:
        """
        sample_stats.py \\
            --sample_name ${sampleName} \\
            --flagstat_file ${flagStats}  \\
            --samtoolsstats_file ${samtoolsStats} \\
            --wgsmetrics_file ${wgsMetrics} \\
            --ntmfraction_file ${ntmFraction} \\
            --cutoff_median_coverage ${params.cutoff_median_coverage} \\
            --cutoff_breadth_of_coverage ${params.cutoff_breadth_of_coverage} \\
            --cutoff_ntm_fraction ${params.cutoff_ntm_fraction}
        """

}
