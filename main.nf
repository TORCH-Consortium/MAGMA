nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================

include { MAP_WF } from './workflows/map.nf'
include { QUANTTB_QUANT } from './modules/quanttb/quant.nf' addParams( params.QUANTTB_QUANT )
include { LOFREQ_CALL_NTM } from './modules/lofreq/call_ntm.nf' addParams( params.LOFREQ_CALL_NTM )

//================================================================================
// Prepare channels
//================================================================================


//   0     1       2       3    4  5     6      7       8
// Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence


reads_ch = Channel.fromPath("${projectDir}/data/mock_data/input_samplesheet.csv")
        .splitCsv(header: false, skip: 1)
        .map { row -> {
                    study = row[0]
                    sample = row[1]
                    library = row[2]
                    attempt = row[3]
                    read1 = row[4]
                    read2 = row[5]
                    flowcell = row[6]
                    lane = row[7]
                    index_sequence = row[8]

            unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"

            //NOTE: Accomodate single/multi reads
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

    // MAP_WF(reads_ch)
    // QUANTTB_QUANT(reads_ch)

    lofreq_all_ntm_ch = Channel.of(["SCIENSANO.07-0789.L1.A1.1.1.1", "${projectDir}/results_mock/bwa/mem/SCIENSANO.07-0789.L1.A1.1.1.1.sorted_reads.bam"])

    LOFREQ_CALL_NTM(lofreq_all_ntm_ch, params.ref_fasta)
}



//================================================================================
// Main workflow
//================================================================================

workflow {

    //DONE
    MAP_WF(reads_ch)
    QUANTTB(reads_ch)

    //TODO
    CALL_WF(MAP_WF.out)

    MERGE_WF(QUANTTB.out, CALL_WF.out)

}
