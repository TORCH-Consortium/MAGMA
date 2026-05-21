process UTILS_FORMAT_DR_REGION_COVERAGE {
    tag "${sampleName}"
    label 'cpu_2_memory_2'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bedcovTsv)

    output:
        tuple val(sampleName), path("${sampleName}.dr_region_coverage.tsv"), emit: coverage

    script:
        """
        coverage_stats_dr_regions.py \\
            --sample-name ${sampleName} \\
            --bedcov ${bedcovTsv} \\
            --gene-name-column 5 \\
            --output ${sampleName}.dr_region_coverage.tsv
        """
}
