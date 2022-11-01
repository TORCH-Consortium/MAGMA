include { BWA_MEM } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)

workflow MAP_WF {
    take:
        approved_samples_ch

    main:

    //TODO: Rename and remove the approved
        BWA_MEM__APPROVED(approved_samples_ch,
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

}
