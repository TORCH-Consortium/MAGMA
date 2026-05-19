process UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(sampleStats), path(drCoverage)

    output:
        path("${sampleName}.stats.tsv"), emit: stats

    script:
        """
        merge_sample_stats_dr_coverage.py \\
            --sample-name ${sampleName} \\
            --sample-stats ${sampleStats} \\
            --dr-coverage ${drCoverage} \\
            --output ${sampleName}.stats.tsv
        """
}
