include { FASTQC } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { QUANTTB_QUANT } from '../modules/quanttb/quant.nf' addParams( params.QUANTTB_QUANT )
include { UTILS_QUANTTB_SAMPLE_QC } from '../modules/utils/quanttb_sample_qc.nf' addParams( params.UTILS_QUANTTB_SAMPLE_QC )
include { UTILS_QUANTTB_COHORT_STATS } from '../modules/utils/quanttb_cohort_stats.nf' addParams( params.UTILS_QUANTTB_COHORT_STATS )


workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)

        //TODO: Add FASTQ validator

        //FIXME: Accommodate Lennert's script 
        //QUANTTB_QUANT(reads_ch)

        /* UTILS_QUANTTB_SAMPLE_QC(QUANTTB_QUANT.out.quanttb_report_tuple, */
        /*                         params.rel_abundance_cutoff, */
        /*                         QUANTTB_QUANT.out.samplereads_tuple) */


        /* UTILS_QUANTTB_COHORT_STATS( */
        /*     UTILS_QUANTTB_SAMPLE_QC.out.quanttb_sample_qc_ch.collect() */
        /* ) */

    emit:
        approved_samples_ch = UTILS_QUANTTB_COHORT_STATS.out.approved_samples_tsv
                                                        .splitCsv(header: false, skip: 1, sep: '\t')
                                                        .map { row -> {
                                                                    def derived_sample_name =	row[-1]
                                                                    return tuple("${derived_sample_name}")
                                                                }
                                                            }
                                                        .join(reads_ch)

        reports_fastqc_ch =  FASTQC.out.collect()

}
