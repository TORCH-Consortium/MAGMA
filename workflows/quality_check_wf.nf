include { FASTQC } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { QUANTTB_QUANT } from '../modules/quanttb/quant.nf' addParams( params.QUANTTB_QUANT )
include { UTILS_QUANTTB_SAMPLE_QC } from '../modules/utils/quanttb_sample_qc.nf' addParams( params.UTILS_QUANTTB_SAMPLE_QC )
include { UTILS_QUANTTB_COHORT_STATS } from '../modules/utils/quanttb_cohort_stats.nf' addParams( params.UTILS_QUANTTB_COHORT_STATS )


workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)

        QUANTTB_QUANT(reads_ch)

        UTILS_QUANTTB_SAMPLE_QC(QUANTTB_QUANT.out.quanttb_report_tuple,
                                params.rel_abundance_cutoff,
                                QUANTTB_QUANT.out.samplereads_tuple)


        UTILS_QUANTTB_COHORT_STATS(
            UTILS_QUANTTB_SAMPLE_QC.out.quanttb_sample_qc_ch.collect()
        )

    emit:
        approved_samples_ch = UTILS_QUANTTB_COHORT_STATS.out.approved_samples

        rejected_samples_ch = UTILS_QUANTTB_COHORT_STATS.out.rejected_samples

        reports_fastqc_ch =  FASTQC.out.collect()

}
