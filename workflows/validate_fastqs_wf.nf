include { FASTQ_VALIDATOR } from '../modules/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )
include { UTILS_FASTQ_COHORT_VALIDATION } from '../modules/utils/fastq_cohort_validation.nf' addParams ( params.UTILS_FASTQ_COHORT_VALIDATION  )
include { SAMPLESHEET_VALIDATION } from '../modules/utils/samplesheet_validation.nf' 

workflow VALIDATE_FASTQS_WF {
    take:
         samplesheet

    main:
        SAMPLESHEET_VALIDATION(samplesheet)

        //NOTE: Expected structure of input CSV samplesheet
        //   0     1       2       3    4  5     6      7       8
        // Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence

        reads_ch = SAMPLESHEET_VALIDATION.out
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

                    //NOTE: Platform is hard-coded to illumina
                    bam_rg_string ="@RG\\tID:${flowcell}.${lane}\\tSM:${study}.${sample}\\tPL:illumina\\tLB:lib${library}\\tPU:${flowcell}.${lane}.${index_sequence}"

                    unique_sample_id = "${study}.${sample}.L${library}.A${attempt}.${flowcell}.${lane}.${index_sequence}"


                    println(row)

                    //Accomodate single/multi reads
                    if (read1 && read2) {

                        return [unique_sample_id, bam_rg_string, [file(read1, checkIfExists: true), file(read2, checkIfExists: true)]]

                    } else {

                        return [unique_sample_id, bam_rg_string, [file(read1, checkIfExists: true)]]

                    }
                }
            }

        reads_ch.view()

        FASTQ_VALIDATOR(reads_ch)

        UTILS_FASTQ_COHORT_VALIDATION( FASTQ_VALIDATOR.out.check_result.collect() )

    emit:

        passed_fastqs_ch = UTILS_FASTQ_COHORT_VALIDATION.out.passed_fastqs
                                                        .splitCsv(header: false, sep: '\t')
                                                        .map { row -> { row[0] } }
                                                        .join(reads_ch)


}
