include { QUANTTB_QUANT } from '../modules/quanttb/quant.nf' addParams( params.QUANTTB_QUANT )
include { UTILS_QUANTTB_SAMPLE_QC } from '../modules/utils/quanttb_sample_qc.nf' addParams( params.UTILS_QUANTTB_SAMPLE_QC )
include { UTILS_QUANTTB_COHORT_STATS } from '../modules/utils/quanttb_cohort_stats.nf' addParams( params.UTILS_QUANTTB_COHORT_STATS )


workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        QUANTTB_QUANT(reads_ch)

        UTILS_QUANTTB_SAMPLE_QC(QUANTTB_QUANT.out.quanttb_report_tuple,
                                params.rel_abundance_cutoff,
                                QUANTTB_QUANT.out.samplereads_tuple)


        UTILS_QUANTTB_COHORT_STATS(
            UTILS_QUANTTB_SAMPLE_QC.out.quanttb_sample_qc_ch.collect()
        )

    emit:
        approved_samples = UTILS_QUANTTB_SAMPLE_QC.out.qc_samplereads_tuple
                            .filter {
                                //TODO
                                //find the location of the qc file in tuple
                                //use file.text to read it's text and split on ","
                                //see if the qc threshold test was passed
                                // return true


                                //     if(relabundance_threshold_met == "1") {
                                //         return tuple("${derived_sample_name}")
                                //     }

                            }
                            .map {
                                //TODO
                                //reshape the channel shape here
                                //drop the QC file
                            }
                            .collect()
                            .view()

}
