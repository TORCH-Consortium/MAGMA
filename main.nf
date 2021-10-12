nextflow.enable.dsl = 2


//================================================================================
// Include sub-workflows/modules and (soft) override workflow-level parameters
//================================================================================


//================================================================================
// Prepare channels
//================================================================================

// unique_sample_id = '{study}.{sample}.L{library}.A{attempt}.{flowcell}.{lane}.{index_sequence}'


//================================================================================
// TEST workflow
//================================================================================
workflow TEST {

//  0      1       2       3    4  5     6      7      8
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
            tuple(unique_sample_id, file(read1), file(read2))
            }
        }


reads_ch.view()
}



//================================================================================
// Main workflow
//================================================================================

workflow {

    MAP_WF()

    QUANTTB(MAP_WF.out)
    CALL_WF(MAP_WF.out)

    MERGE_WF(QUANTTB.out, CALL_WF.out)

}
