include { FASTQC } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { MULTIQC } from '../modules/multiqc/multiqc.nf' addParams (params.MULTIQC)
include { BWA_MEM } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)

workflow MAP_WF {
    take:
        approved_reads_ch

    main:

        FASTQC(reads_ch)

        MULTIQC(FASTQC.out.collect())

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
