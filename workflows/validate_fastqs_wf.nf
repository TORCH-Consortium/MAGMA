include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )
include { UTILS_FASTQ_COHORT_VALIDATION } from '../modules/utils/fastq_cohort_validation.nf' addParams ( params.UTILS_FASTQ_COHORT_VALIDATION  )

workflow VALIDATE_FASTQS_WF {
    take:
         reads_ch

    main:

        //FIXME: Add the samplesheet validator process for samplesheet_validation.py script

        FASTQ_VALIDATOR(reads_ch)

        UTILS_FASTQ_COHORT_VALIDATION( FASTQ_VALIDATOR.out.check_result.collect() )

    emit:

        passed_fastqs_ch = UTILS_FASTQ_COHORT_VALIDATION.out.passed_fastqs
                                                        .splitCsv(header: false, sep: '\t')
                                                        .map { row -> { row[0] } }
                                                        .join(reads_ch)


}
