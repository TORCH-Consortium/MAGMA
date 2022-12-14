nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================

include { CALL_WF } from './workflows/call_wf.nf'
include { VALIDATE_FASTQS_WF } from './workflows/validate_fastqs_wf.nf'
include { MULTIPLE_INFECTIONS_WF } from './workflows/multiple_infections_wf.nf'
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

        validated_reads_ch = VALIDATE_FASTQS_WF(params.input_samplesheet)

        QUALITY_CHECK_WF(validated_reads_ch)

        MAP_WF( QUALITY_CHECK_WF.out.approved_samples_ch,
                QUALITY_CHECK_WF.out.rejected_samples_ch )

        MULTIPLE_INFECTIONS_WF(MAP_WF.out.rejected_sorted_reads_ch)

        CALL_WF(MAP_WF.out.approved_sorted_reads_ch)

        collated_gvcfs_ch = CALL_WF.out.gvcf_ch
            .flatten()
            .collate(3)
            .dump(tag:'XBS-NF-LOG collated_gvcfs_ch : ', pretty: true)

        sample_stats_ch = CALL_WF.out.cohort_stats_tsv
            .splitCsv(header: false, skip: 1, sep: '\t' )
            .map { row -> [
                    row.first(),           // SAMPLE
                    row.last().toInteger() // ALL_THRESHOLDS_MET
            ]
        }
            .filter { it[1] == 1} // Filter out samples which meet all the thresholds
            .map { [ it[0] ] }
            .dump(tag:'XBS-NF-LOG sample_stats_ch : ', pretty: true)

        selected_gvcfs_ch = collated_gvcfs_ch.join(sample_stats_ch)
            .flatten()
            .filter { it -> { 
                            (it.class.name  == sun.nio.fs.UnixPath) 
                            || (it.class.name == "nextflow.cloud.azure.nio.AzPath") 
                            || (it.class.name == "com.upplication.s3fs.S3Path") 
                            || (it.class.name == "com.google.cloud.storage.contrib.nio.CloudStoragePath") 
                    } }
                .collect()
            .dump(tag:'XBS-NF-LOG selected_gvcfs_ch.flatten().filter().collect() : ', pretty: true)

        MERGE_WF(selected_gvcfs_ch, CALL_WF.out.lofreq_vcf_ch)


        REPORTS_WF(QUALITY_CHECK_WF.out.reports_fastqc_ch)
    }

}

