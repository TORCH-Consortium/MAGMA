include { FASTQC } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)

workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)

    emit:
        reports_fastqc_ch =  FASTQC.out.collect()

}
