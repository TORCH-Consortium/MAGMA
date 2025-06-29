/*
 * Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
 *
 * This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
 *
 * For quick overview of GPL-3 license, please refer
 * https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
 *
 * - You MUST keep this license with original authors in your copy
 * - You MUST acknowledge the original source of this software
 * - You MUST state significant changes made to the original software
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program . If not, see <http://www.gnu.org/licenses/>.
 */
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
