/*
 * Copyright ...
 */

include { UTILS_SAMPLE_STATS as UTILS_RAW_SAMPLE_STATS } from '../../modules/local/utils/sample_stats' addParams( params.UTILS_SAMPLE_STATS )
include { SAMTOOLS_COVERAGE_STATS_DR_REGIONS } from '../../modules/local/samtools/coverage_stats_dr_regions' addParams( params.SAMTOOLS_COVERAGE_STATS_DR_REGIONS )
include { UTILS_FORMAT_DR_REGION_COVERAGE } from '../../modules/local/utils/format_dr_region_coverage' addParams( params.UTILS_FORMAT_DR_REGION_COVERAGE )
include { UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE } from '../../modules/local/utils/merge_sample_stats_dr_coverage' addParams( params.UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE )

workflow UTILS_SAMPLE_STATS {

    take:
        ch_sample_stats_input
        ch_bam

    main:
        ch_dr_regions = Channel.value(
            file("${projectDir}/resources/regions/tbprofiler_whov2plus_genes.bed")
        )

        UTILS_RAW_SAMPLE_STATS(ch_sample_stats_input)

        SAMTOOLS_COVERAGE_STATS_DR_REGIONS(
            ch_bam,
            ch_dr_regions
        )
        
        UTILS_FORMAT_DR_REGION_COVERAGE(
            SAMTOOLS_COVERAGE_STATS_DR_REGIONS.out.bedcov
        )
        
        ch_merged_input = UTILS_RAW_SAMPLE_STATS.out
            .join(UTILS_FORMAT_DR_REGION_COVERAGE.out.coverage)
            .map { sampleName, sampleStats, drCoverage ->
                tuple(sampleName, sampleStats, drCoverage)
            }

UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE(ch_merged_input)

    emit:
        UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE.out
}
