nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================

include { MAP_WF } from './workflows/map_wf.nf'
include { QUANTTB_QUANT } from './modules/quanttb/quant.nf' addParams( params.QUANTTB_QUANT )
include { CALL_WF } from './workflows/call_wf.nf'

//================================================================================
// Prepare channels
//================================================================================


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

            //TODO: Confirm whether '.' can be replaced with '_'
            unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"

            //NOTE: Accomodate single/multi reads
            //TODO: Consider using read1 and read2 identifiers - confirm before refactor.
            if (row[4] && row[5]) {

                // Both read1 and read2 are present
                return tuple(unique_sample_id, tuple(file(read1), file(read2)))

            } else if (row[4]) {

                // Only read1 is present
                return tuple(unique_sample_id, tuple(file(read1)))

            } else {

                // Only read2 is present
                return tuple(unique_sample_id, tuple(file(read2)))

            }
        }
    }




//================================================================================
// TEST workflow
//================================================================================

workflow TEST {
    QUANTTB_QUANT(reads_ch)
    MAP_WF(reads_ch)
    CALL_WF(MAP_WF.out.sorted_reads, QUANTTB_QUANT.out)



    CALL_WF.out.cohort_stats_tsv
        .splitCsv(header: true )
        .view()

    // .map { row -> {
    //                 study           = row[0]
    //                 sample          = row[1]
    //                 library         = row[2]
    //                 attempt         = row[3]
    //                 read1           = row[4]
    //                 read2           = row[5]
    //                 flowcell        = row[6]
    //                 lane            = row[7]
    //                 index_sequence  = row[8]


    //             SAMPLE
    //             AVG_INSERT_SIZE
    //             MAPPED_PERCENTAGE
    //             RAW_TOTAL_SEQS
    //             AVERAGE_QUALITY
    //             QUANTTB_RELATIVE_ABUNDANCE
    //             QUANTTB_DEPTH
    //             MEAN_COVERAGE
    //             SD_COVERAGE
    //             MEDIAN_COVERAGE
    //             MAD_COVERAGE
    //             PCT_EXC_ADAPTER
    //             PCT_EXC_MAPQ
    //             PCT_EXC_DUPE
    //             PCT_EXC_UNPAIRED
    //             PCT_EXC_BASEQ
    //             PCT_EXC_OVERLAP
    //             PCT_EXC_CAPPED
    //             PCT_EXC_TOTAL
    //             PCT_1X
    //             PCT_5X
    //             PCT_10X
    //             PCT_30X
    //             PCT_50X
    //             PCT_100X
    //             NTM_FRACTION
    //             NTM_FRACTION_THRESHOLD_MET
    //             RELATIVE_ABUNDANCE_THRESHOLD_MET
    //             COVERAGE_THRESHOLD_MET
    //             BREADTH_OF_COVERAGE_THRESHOLD_MET
    //             ALL_THRESHOLDS_MET

    //         //TODO: Confirm whether '.' can be replaced with '_'
    //         unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"
    //     }
    // }



    // CALL_WF.out.gvcf_ch.view()

}



//================================================================================
// Main workflow
//================================================================================

workflow {

    //DONE
    QUANTTB_QUANT(reads_ch)
    MAP_WF(reads_ch)
    CALL_WF(MAP_WF.out.sorted_reads)

    //TODO
    // MERGE_WF(QUANTTB.out, CALL_WF.out)

}
