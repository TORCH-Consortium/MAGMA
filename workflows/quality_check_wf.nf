include { FASTQC              } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)

workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)


        def ntmprofilerDb =  params.ntmprofilerDb != "DEFAULT" ?  params.ntmprofilerDb : []



    emit:
        reports_fastqc_ch =  FASTQC.out.collect()

}
