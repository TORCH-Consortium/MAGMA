nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================


include { CALL_WF } from './workflows/call_wf.nf'
include { VALIDATE_FASTQS_WF } from './workflows/validate_fastqs_wf.nf'
include { MAP_WF } from './workflows/map_wf.nf'
include { MERGE_WF } from './workflows/merge_wf.nf'
include { MINOR_VARIANTS_ANALYSIS_WF } from './workflows/minor_variants_analysis_wf.nf'
include { QUALITY_CHECK_WF } from './workflows/quality_check_wf.nf'
include { REPORTS_WF } from './workflows/reports_wf.nf'
include { STRUCTURAL_VARIANTS_ANALYSIS_WF } from './workflows/structural_variants_analysis_wf.nf'
include { UTILS_MERGE_COHORT_STATS } from "./modules/utils/merge_cohort_stats.nf" addParams ( params.UTILS_MERGE_COHORT_STATS )

//================================================================================
// Main workflow
//================================================================================

workflow {

    if (params.only_validate_fastqs) {

        VALIDATE_FASTQS_WF(params.input_samplesheet)

    } else if (params.only_validate_and_qc) {

        validated_reads_ch = VALIDATE_FASTQS_WF( params.input_samplesheet )

        QUALITY_CHECK_WF( validated_reads_ch )

    } else if (params.skip_merge_analysis) {

        validated_reads_ch = VALIDATE_FASTQS_WF( params.input_samplesheet )

        QUALITY_CHECK_WF( validated_reads_ch )

        MAP_WF( validated_reads_ch )

        CALL_WF( MAP_WF.out.sorted_reads_ch )

        MINOR_VARIANTS_ANALYSIS_WF(CALL_WF.out.reformatted_lofreq_vcfs_tuple_ch)

        UTILS_MERGE_COHORT_STATS ( MINOR_VARIANTS_ANALYSIS_WF.out.approved_samples_ch,
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

        STRUCTURAL_VARIANTS_ANALYSIS_WF ( validated_reads_ch, all_samples_ch )


    //FIXME IMPLEMENT ANOTHER LOGIC FOR --only_merge_analysis
    /* else if { */
    /*         //--cohort_data  --only_merge_wf */
    /* } */

    } else {

        validated_reads_ch = VALIDATE_FASTQS_WF( params.input_samplesheet )

        QUALITY_CHECK_WF( validated_reads_ch )

        MAP_WF( validated_reads_ch )

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

        STRUCTURAL_VARIANTS_ANALYSIS_WF ( validated_reads_ch, all_samples_ch )


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
                    MINOR_VARIANTS_ANALYSIS_WF.out.minor_variants_results_ch,
                    MERGE_WF.out.major_variants_results_ch,
                    STRUCTURAL_VARIANTS_ANALYSIS_WF.out.structural_variants_results_ch )

    }

}

