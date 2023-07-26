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
            --median_coverage_cutoff ${params.median_coverage_cutoff} \\
            --breadth_of_coverage_cutoff ${params.breadth_of_coverage_cutoff} \\
            --rel_abundance_cutoff ${params.rel_abundance_cutoff} \\
            --ntm_fraction_cutoff ${params.ntm_fraction_cutoff}
        """

}
