include { BWA_MEM as BWA_MEM__APPROVED } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)
include { BWA_MEM as BWA_MEM__REJECTED } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM__REJECTED)

workflow MAP_WF {
    take:
        approved_samples_ch
        rejected_samples_ch

    main:

        BWA_MEM__APPROVED(approved_samples_ch,
                          params.ref_fasta,
                          [params.ref_fasta_dict,
                           params.ref_fasta_amb,
                           params.ref_fasta_ann,
                           params.ref_fasta_bwt,
                           params.ref_fasta_fai,
                           params.ref_fasta_pac,
                           params.ref_fasta_sa])


        BWA_MEM__REJECTED(rejected_samples_ch,
                          params.ref_fasta,
                          [params.ref_fasta_dict,
                           params.ref_fasta_amb,
                           params.ref_fasta_ann,
                           params.ref_fasta_bwt,
                           params.ref_fasta_fai,
                           params.ref_fasta_pac,
                           params.ref_fasta_sa])


    emit:
        approved_sorted_reads_ch = BWA_MEM__APPROVED.out
        rejected_sorted_reads_ch = BWA_MEM__REJECTED.out

}
