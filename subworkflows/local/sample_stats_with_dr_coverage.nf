/*
 * Copyright ...
 */

include { UTILS_SAMPLE_STATS as UTILS_RAW_SAMPLE_STATS } from '../../modules/local/utils/sample_stats'
include { SAMTOOLS_COVERAGE_STATS_DR_REGIONS } from '../../modules/local/samtools/coverage_stats_dr_regions'
include { UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE } from '../../modules/local/utils/merge_sample_stats_dr_coverage'

workflow UTILS_SAMPLE_STATS {

    take:
        ch_sample_stats_input
        ch_bam

    main:
        ch_dr_regions = Channel.value(
            file("${projectDir}/resources/regions/WHO_Tier1_Tier2_DR.list")
        )

        UTILS_RAW_SAMPLE_STATS(ch_sample_stats_input)

        SAMTOOLS_COVERAGE_STATS_DR_REGIONS(
            ch_bam,
            ch_dr_regions
        )

        ch_merged_input = UTILS_RAW_SAMPLE_STATS.out
            .join(SAMTOOLS_COVERAGE_STATS_DR_REGIONS.out)
            .map { sampleName, sampleStats, drCoverage ->
                tuple(sampleName, sampleStats, drCoverage)
            }
        
        UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE(ch_merged_input)

    emit:
        UTILS_MERGE_SAMPLE_STATS_DR_COVERAGE.out
}
