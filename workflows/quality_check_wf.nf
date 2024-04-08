include { FASTQC              } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { NTMPROFILER_PROFILE } from '../modules/ntmprofiler/profile.nf' addParams (params.NTMPROFILER_PROFILE)
include { NTMPROFILER_COLLATE } from '../modules/ntmprofiler/collate.nf' addParams (params.NTMPROFILER_COLLATE)

workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)

        NTMPROFILER_PROFILE( reads_ch )

        NTMPROFILER_COLLATE( params.vcf_name,
                             NTMPROFILER_PROFILE.out.resistance_json.collect() )


    emit:
        reports_fastqc_ch =  FASTQC.out.collect()

}
