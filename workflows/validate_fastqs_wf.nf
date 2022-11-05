include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )

workflow VALIDATE_FASTQS {
    take:
         reads_ch

    main:
        FASTQ_VALIDATOR(reports_fastqc_ch)

}
