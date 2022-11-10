nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================

include { CALL_WF } from './workflows/call_wf.nf'
include { VALIDATE_FASTQS_WF } from './workflows/validate_fastqs_wf.nf'
include { MAP_WF } from './workflows/map_wf.nf'
include { MERGE_WF } from './workflows/merge_wf.nf'
include { QUALITY_CHECK_WF } from './workflows/quality_check_wf.nf'
include { REPORTS_WF } from './workflows/reports_wf.nf'


//================================================================================
// Main workflow
//================================================================================

workflow {

    if (params.only_validate_fastqs) {

        VALIDATE_FASTQS_WF(params.input_samplesheet)

    } else {

        validated_reads_ch = VALIDATE_FASTQS_WF( params.input_samplesheet )

        QUALITY_CHECK_WF( validated_reads_ch )

        MAP_WF( validated_reads_ch )

        CALL_WF( MAP_WF.out.sorted_reads_ch )

        //FIXME
        MINOR_VARIANT_ANALYSIS_WF (CALL_WF.FIXME)

        //FIXME Take the results of MINOR_VARIANT_ANALYSIS_WF analysis for the samples which are approved 
        //and combine with existing filtering process.
        //Combine the results with those of cohort stats and then do the filtering

        collated_gvcfs_ch = CALL_WF.out.gvcf_ch
            .flatten()
            .collate(3)
            // .view(it -> "\n\n XBS-NF-LOG collated_gvcfs_ch : $it \n\n")

        sample_stats_ch = CALL_WF.out.cohort_stats_tsv
            .splitCsv(header: false, skip: 1, sep: '\t' )
            .map { row -> [
                    row.first(),           // SAMPLE
                    row.last().toInteger() // ALL_THRESHOLDS_MET
            ]
        }
            .filter { it[1] == 1} // Filter out samples which meet all the thresholds
            .map { [ it[0] ] }
            // .view("\n\n XBS-NF-LOG sample_stats_ch : $it \n\n")

        selected_gvcfs_ch = collated_gvcfs_ch.join(sample_stats_ch)
            .flatten()
            .filter { it.class  == sun.nio.fs.UnixPath }
            // .view("\n\n XBS-NF-LOG selected_gvcfs_ch : $it \n\n")


        MERGE_WF(selected_gvcfs_ch.collect(), CALL_WF.out.reformatted_lofreq_vcf_ch)


        REPORTS_WF(QUALITY_CHECK_WF.out.reports_fastqc_ch)
    }

}

