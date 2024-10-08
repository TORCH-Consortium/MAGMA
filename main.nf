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
nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================


include { CALL_WF } from './workflows/call_wf.nf'
include { VALIDATE_FASTQS_WF } from './workflows/validate_fastqs_wf.nf'
include { MAP_WF } from './workflows/map_wf.nf'
include { MERGE_WF } from './workflows/merge_wf.nf'
include { MINOR_VARIANTS_ANALYSIS_WF } from './workflows/minor_variants_analysis_wf.nf'
//include { MULTIQC AS MULTIQC_FASTQS } from '../modules/multiqc/multiqc.nf' addParams (params.MULTIQC_FASTQS)
include { QUALITY_CHECK_WF } from './workflows/quality_check_wf.nf'
include { REPORTS_WF } from './workflows/reports_wf.nf'
include { SAMPLESHEET_VALIDATION } from './modules/utils/samplesheet_validation.nf'  addParams ( params.SAMPLESHEET_VALIDATION )
include { STRUCTURAL_VARIANTS_ANALYSIS_WF } from './workflows/structural_variants_analysis_wf.nf'
include { UTILS_MERGE_COHORT_STATS } from "./modules/utils/merge_cohort_stats.nf" addParams ( params.UTILS_MERGE_COHORT_STATS )

//================================================================================
// Main workflow
//================================================================================

workflow {

    if (params.only_validate_fastqs) {

        SAMPLESHEET_VALIDATION( params.input_samplesheet )

        VALIDATE_FASTQS_WF( SAMPLESHEET_VALIDATION.out.validated_samplesheet , SAMPLESHEET_VALIDATION.out.status )

        QUALITY_CHECK_WF( VALIDATE_FASTQS_WF.out.approved_fastqs_ch )

        //TODO: Add modules for generating fastq stats and then capturing them in the MultiQC image
        //MULTIQC_FASTQS( QUALITY_CHECK_WF.out.reports_fastqc_ch )

    } else  {

        SAMPLESHEET_VALIDATION(params.input_samplesheet)


        VALIDATE_FASTQS_WF( SAMPLESHEET_VALIDATION.out.validated_samplesheet , SAMPLESHEET_VALIDATION.out.status )

        QUALITY_CHECK_WF( VALIDATE_FASTQS_WF.out.approved_fastqs_ch )


        MAP_WF( VALIDATE_FASTQS_WF.out.approved_fastqs_ch  )

        CALL_WF( MAP_WF.out.sorted_reads_ch )

        //NOTE: Samples implicitly get filtered in BCFTOOLS_MERGE if they don't have any identified variants
        MINOR_VARIANTS_ANALYSIS_WF(CALL_WF.out.reformatted_lofreq_vcfs_tuple_ch)

        UTILS_MERGE_COHORT_STATS( MINOR_VARIANTS_ANALYSIS_WF.out.approved_samples_ch,
                                  MINOR_VARIANTS_ANALYSIS_WF.out.rejected_samples_ch,
                                  CALL_WF.out.cohort_stats_tsv )


        all_samples_ch = UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch
                                .splitCsv(header: false, skip: 1, sep: '\t' )
                                .map { row -> [
                                        row.first(),           // SAMPLE
                                        row.last().toInteger() // ALL_THRESHOLDS_MET
                                        ]
                                    }
                                .map { [ it[0] ] }
                                //.dump(tag:'MERGE_WF: all_samples_ch', pretty: true)

        STRUCTURAL_VARIANTS_ANALYSIS_WF ( VALIDATE_FASTQS_WF.out.approved_fastqs_ch, all_samples_ch )


        if (!params.skip_merge_analysis) {

            approved_samples_ch = UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch
                                    .splitCsv(header: false, skip: 1, sep: '\t' )
                                    .map { row -> [
                                            row.first(),           // SAMPLE
                                            row.last().toInteger() // ALL_THRESHOLDS_MET
                                            ]
                                        }
                                    .filter { it[1] == 1} // Filter out samples which meet all the thresholds
                                    .map { [ it[0] ] }
                                    //.dump(tag:'MERGE_WF: approved_samples_ch', pretty: true)


            MERGE_WF( CALL_WF.out.gvcf_ch,
                      CALL_WF.out.reformatted_lofreq_vcfs_tuple_ch,
                      approved_samples_ch )


            REPORTS_WF( QUALITY_CHECK_WF.out.reports_fastqc_ch,
                        UTILS_MERGE_COHORT_STATS.out.merged_cohort_stats_ch,
                        MERGE_WF.out.major_variants_results_ch,
                        MINOR_VARIANTS_ANALYSIS_WF.out.minor_variants_results_ch,
                        STRUCTURAL_VARIANTS_ANALYSIS_WF.out.structural_variants_results_ch )

        }
    }
}
