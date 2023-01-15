include { BWA_MEM } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)

workflow MAP_WF {
    take:
        samples_ch

    main:

        BWA_MEM(samples_ch,
                  params.ref_fasta,
                  [params.ref_fasta_dict,
                   params.ref_fasta_amb,
                   params.ref_fasta_ann,
                   params.ref_fasta_bwt,
                   params.ref_fasta_fai,
                   params.ref_fasta_pac,
                   params.ref_fasta_sa])


    emit:
        sorted_reads_ch = BWA_MEM.out

}
