//FIXME Replace the name of FASTQ_VALIDATOR with FASTQ_REPORT
//include { UTILS_FASTQ_REPORT } from '../modules/utils/fastq_report.nf' addParams (params.FASTQ_REPORT)
include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )
include { UTILS_FASTQ_COHORT_VALIDATION } from '../modules/utils/fastq_cohort_validation.nf' addParams ( params.UTILS_FASTQ_COHORT_VALIDATION  )

workflow VALIDATE_FASTQS_WF {
    take:
         samplesheet
         ready

    main:

    //NOTE: Expected structure of input CSV samplesheet
    //   0     1       2       3    4  5     6      7       8
    // Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence



        fastqs_ch = samplesheet
                    .splitCsv(header: false, skip: 1)
                    .map { row -> {
                                study           = row[0]
                                sample          = row[1]
                                library         = row[2]
                                attempt         = row[3]
                                read1           = row[4]
                                read2           = row[5]
                                flowcell        = row[6]
                                lane            = row[7]
                                index_sequence  = row[8]
                                magma_sample_name = row[9]
                                magma_bam_rg_string = row[10]


                //Accomodate single/multi reads
                if (read1 && read2) {

                    return [[id: magma_sample_name, paired: true, bam_rg_string:magma_bam_rg_string ], [file(read1, checkIfExists: true), file(read2, checkIfExists: true)]]

                } else {

                    return [[id: magma_sample_name, paired: true, bam_rg_string:magma_bam_rg_string ],  [file(read1, checkIfExists: true)]]

                    }
                }
            }.transpose()


    FASTQ_VALIDATOR( fastqs_ch, ready )


    UTILS_FASTQ_COHORT_VALIDATION( FASTQ_VALIDATOR.out.fastq_report.collect(), samplesheet )



    emit:

        passed_fastqs_ch = UTILS_FASTQ_COHORT_VALIDATION.out.passed_fastqs
                                                            .splitText()
                                                            .join(reads_ch)
                                                            .view()
}
