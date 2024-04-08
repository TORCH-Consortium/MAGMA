include { FASTQC              } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { NTMPROFILER_PROFILE } from '../modules/ntmprofiler/profile.nf' addParams (params.NTMPROFILER_PROFILE)

workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)

        NTMPROFILER_PROFILE( reads_ch )


    emit:
        reports_fastqc_ch =  FASTQC.out.collect()

}
