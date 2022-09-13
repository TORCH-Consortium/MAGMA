process UTILS_COHORT_STATS {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("sample_stats/*")

    output:
        path("*.cohort_stats.tsv")

    shell:
        '''
        echo -e "SAMPLE\tAVG_INSERT_SIZE\tMAPPED_PERCENTAGE\tRAW_TOTAL_SEQS\tAVERAGE_BASE_QUALITY\tMEAN_COVERAGE\tSD_COVERAGE\tMEDIAN_COVERAGE\tMAD_COVERAGE\tPCT_EXC_ADAPTER\tPCT_EXC_MAPQ\tPCT_EXC_DUPE\tPCT_EXC_UNPAIRED\tPCT_EXC_BASEQ\tPCT_EXC_OVERLAP\tPCT_EXC_CAPPED\tPCT_EXC_TOTAL\tPCT_1X\tPCT_5X\tPCT_10X\tPCT_30X\tPCT_50X\tPCT_100X\tMAPPED_NTM_FRACTION_16S\tMAPPED_NTM_FRACTION_16S_THRESHOLD_MET\tCOVERAGE_THRESHOLD_MET\tBREADTH_OF_COVERAGE_THRESHOLD_MET\tALL_THRESHOLDS_MET" > !{params.vcf_name}.cohort_stats.tsv
        cat sample_stats/*tsv >> !{params.vcf_name}.cohort_stats.tsv
        '''
}
