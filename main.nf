nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================

include { QUALITY_CHECK_WF } from './workflows/quality_check_wf.nf'
include { MAP_WF } from './workflows/map_wf.nf'
include { CALL_WF } from './workflows/call_wf.nf'
include { MERGE_WF } from './workflows/merge_wf.nf'

//================================================================================
// Prepare channels
//================================================================================


//NOTE: Expected structure of input CSV samplesheet
//   0     1       2       3    4  5     6      7       8
// Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence


reads_ch = Channel.fromPath(params.input_samplesheet)
        .splitCsv(header: false, skip: 1)
        .map { row -> {
                    study           = row[0]
                    sample          = row[1]
                    library         = row[2]
                    attempt         = row[3]
                    read1           = row[4]
                    read2           = row[5]
                    flowcell        = row[6]
                    lane            = row[7]
                    index_sequence  = row[8]

            bam_rg_string ="@RG\\tID:${flowcell}.${lane}\\tSM:${study}.${sample}\\tPL:illumina\\tLB:lib${library}\\tPU:${flowcell}.${lane}.${index_sequence}"

            unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"

            //Accomodate single/multi reads
            if (read1 && read2) {

                return tuple(unique_sample_id, bam_rg_string, tuple(file(read1), file(read2)))

            } else if (read1) {

                return tuple(unique_sample_id, bam_rg_string, tuple(file(read1)))

            } else {

                return tuple(unique_sample_id, bam_rg_string, tuple(file(read2)))

            }
        }
    }




//================================================================================
// TEST workflow
//================================================================================

workflow TEST {

    QUALITY_CHECK_WF(reads_ch)

    // MAP_WF(QUALITY_CHECK_WF.out)

    // CALL_WF(MAP_WF.out.sorted_reads)

    // collated_gvcfs_ch = CALL_WF.out.gvcf_ch.flatten().collate(3)

    // // collated_gvcfs_ch.view()

    // sample_stats_ch = CALL_WF.out.cohort_stats_tsv
    //     .splitCsv(header: false, skip: 1, sep: '\t' )
    //     .map { row -> [
    //             row.first(),           // SAMPLE
    //             row.last().toInteger() // ALL_THRESHOLDS_MET
    //      ]
    // }
    // .filter { it[1] == 1} // Filter out samples which meet all the thresholds
    // .map { [ it[0] ] }
    // // .view()


    // selected_gvcfs_ch = collated_gvcfs_ch.join(sample_stats_ch)
    //     .flatten()
    //     .filter { it.class  == sun.nio.fs.UnixPath }
    //     // .view()


    // MERGE_WF(selected_gvcfs_ch.collect(), CALL_WF.out.lofreq_vcf_ch)

}



//================================================================================
// Main workflow
//================================================================================

workflow {

    QUALITY_CHECK_WF(reads_ch)

    MAP_WF(QUALITY_CHECK_WF.out)

    CALL_WF(MAP_WF.out.sorted_reads)

    collated_gvcfs_ch = CALL_WF.out.gvcf_ch.flatten().collate(3)

    // collated_gvcfs_ch.view()

    sample_stats_ch = CALL_WF.out.cohort_stats_tsv
        .splitCsv(header: false, skip: 1, sep: '\t' )
        .map { row -> [
                row.first(),           // SAMPLE
                row.last().toInteger() // ALL_THRESHOLDS_MET
         ]
    }
    .filter { it[1] == 1} // Filter out samples which meet all the thresholds
    .map { [ it[0] ] }
    // .view()


    selected_gvcfs_ch = collated_gvcfs_ch.join(sample_stats_ch)
        .flatten()
        .filter { it.class  == sun.nio.fs.UnixPath }
        // .view()


    MERGE_WF(selected_gvcfs_ch.collect(), CALL_WF.out.lofreq_vcf_ch)

}
