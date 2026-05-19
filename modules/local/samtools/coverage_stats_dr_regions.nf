process SAMTOOLS_COVERAGE_STATS_DR_REGIONS {
    tag "${sampleName}"

    input:
        tuple val(sampleName), path(bam), path(bai)
        path regions

    output:
        tuple val(sampleName), path("${sampleName}.dr_region_coverage.tsv"), emit: coverage

    script:
        """
        samtools bedcov ${regions} ${bam} > ${sampleName}.dr_region_coverage.long.tsv

        coverage_stats_dr_regions.py \\
            --sample-name ${sampleName} \\
            --bedcov ${sampleName}.dr_region_coverage.long.tsv \\
            --output ${sampleName}.dr_region_coverage.tsv
        """
}
