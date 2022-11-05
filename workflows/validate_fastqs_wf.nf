include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )

workflow VALIDATE_FASTQS_WF {
    take:
         reads_ch

    main:
        FASTQ_VALIDATOR(reads_ch)

}
