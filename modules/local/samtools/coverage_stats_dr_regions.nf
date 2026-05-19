process SAMTOOLS_COVERAGE_STATS_DR_REGIONS {
    tag "${sampleName}"
    label 'cpu_2_memory_2'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bam), path(bai)
        path regions

    output:
        tuple val(sampleName), path("${sampleName}.dr_region_coverage.long.tsv"), emit: bedcov

    script:
        """
        ${params.samtools_path} bedcov ${regions} ${bam} > ${sampleName}.dr_region_coverage.long.tsv
        """
}
