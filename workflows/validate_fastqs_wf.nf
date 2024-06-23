//FIXME Replace the name of FASTQ_VALIDATOR with FASTQ_REPORT
//include { UTILS_FASTQ_REPORT } from '../modules/utils/fastq_report.nf' addParams (params.FASTQ_REPORT)
include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )
include { UTILS_FASTQ_COHORT_VALIDATION } from '../modules/utils/fastq_cohort_validation.nf' addParams ( params.UTILS_FASTQ_COHORT_VALIDATION  )

workflow VALIDATE_FASTQS_WF {
    take:
         samplesheet_json
         ready

    main:

    //NOTE: Expected structure of input CSV samplesheet
    //   0     1       2       3    4  5     6      7       8
    // Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence



    fastqs_ch = samplesheet_json
        .splitJson()
        .map {
            if (it.R2) {
                [it.magma_sample_name, [it.R1, it.R2]]
            } else {

                [it.magma_sample_name, [it.R1]]
            }
        }.transpose()




    FASTQ_VALIDATOR( fastqs_ch, ready )


    UTILS_FASTQ_COHORT_VALIDATION( FASTQ_VALIDATOR.out.fastq_report.collect(), samplesheet_json )


    approved_fastqs_ch = UTILS_FASTQ_COHORT_VALIDATION.out.magma_analysis_json
                            .splitJson()
                            .filter {it.value.fastqs_approved}
                            .map {
                                    if (it.value.R2) {
                                    [it.value.magma_sample_name, it.value.magma_bam_rg_string, [it.value.R1, it.value.R2]]
                                    } else {
                                        [it.value.magma_sample_name, it.value.magma_bam_rg_string, [it.value.R1]]
                                    }
                            }

     emit:

        approved_fastqs_ch = approved_fastqs_ch


}
