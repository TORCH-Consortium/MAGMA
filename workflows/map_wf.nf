include { FASTQC } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { BWA_MEM } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)

workflow MAP_WF {
    take:
        reads_ch

    main:
        FASTQC(reads_ch)


        //TODO: Can be refactored in next iteration to re-use the reads_ch channel
        bew_mem_rg_ch = Channel.fromPath("${projectDir}/data/mock_data/input_samplesheet.csv")
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
            }
              RG="@RG\tID:${flowcell}.${lane}\tSM:${study}.${sample}\tPL:illumina\tLB:lib${library}\tPU:${flowcell}.${lane}.${index_sequence}"
              //TODO: Confirm whether '.' can be replaced with '_'
              unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"
              return tuple(unique_sample_id, RG)
        }


        BWA_MEM(reads_ch.join(bew_mem_rg_ch), params.ref_fasta)

    emit:
        sorted_reads = BWA_MEM.out

}
