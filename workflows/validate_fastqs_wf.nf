include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )

workflow VALIDATE_FASTQS_WF {
    take:
         reads_ch

    main:

        //FIXME: Add the samplesheet validator process for sample_sheet_validation.py script

        FASTQ_VALIDATOR(reads_ch)

    emit:

        validate_reads_ch = reads_ch

}
