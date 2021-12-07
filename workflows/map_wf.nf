include { FASTQC } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { MULTIQC } from '../modules/multiqc/multiqc.nf' addParams (params.MULTIQC)
include { BWA_MEM } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)

workflow MAP_WF {
    take:
        reads_ch

    main:
        FASTQC(reads_ch)

        //TODO: Activate this later, once the renaming is refactored.
        // MULTIQC(FASTQC.out.collect())

        bew_mem_rg_ch = Channel.fromPath(params.input_samplesheet)
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
            }

              RG="@RG\\tID:${flowcell}.${lane}\\tSM:${study}.${sample}\\tPL:illumina\\tLB:lib${library}\\tPU:${flowcell}.${lane}.${index_sequence}"

              unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"

              return tuple(unique_sample_id, RG)
        }


        BWA_MEM(reads_ch.join(bew_mem_rg_ch),
                params.ref_fasta,
                [params.ref_fasta_dict,
                params.ref_fasta_amb,
                params.ref_fasta_ann,
                params.ref_fasta_bwt,
                params.ref_fasta_fai,
                params.ref_fasta_pac,
                params.ref_fasta_sa])



    emit:
        sorted_reads = BWA_MEM.out

}
