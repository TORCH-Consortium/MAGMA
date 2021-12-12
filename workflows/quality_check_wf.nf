include { QUANTTB_QUANT } from '../modules/quanttb/quant.nf' addParams( params.QUANTTB_QUANT )
include { UTILS_QUANTTB_SAMPLE_QC } from '../modules/utils/quanttb_sample_qc.nf' addParams( params.UTILS_QUANTTB_SAMPLE_QC )
include { UTILS_QUANTTB_COHORT_STATS } from '../modules/utils/quanttb_cohort_stats.nf' addParams( params.UTILS_QUANTTB_COHORT_STATS )


workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        QUANTTB_QUANT(reads_ch)

        UTILS_QUANTTB_SAMPLE_QC(QUANTTB_QUANT.out.quanttb_report_tuple,
                                params.rel_abundance_cutoff)

        UTILS_QUANTTB_COHORT_STATS(
            UTILS_QUANTTB_SAMPLE_QC.out.collect()
        )


        // approved_samples_ch = Channel.fromPath("results/quanttb/cohort_stats/joint.quanttb_cohort_stats.tsv")
        approved_samples_ch = UTILS_QUANTTB_COHORT_STATS.out
            .splitCsv(header: false, skip: 1)
            .map{
                row -> {
                    sample_name = row[0]
                    relabundance_threshold_met = row[3]
                    derived_sample_name = row[-1]
                    }

                if(relabundance_threshold_met == "1") {
                    return tuple("${derived_sample_name}")
                }
            }
            .toList()
            .collect()
            .view()


        // temp_reads_ch = reads_ch.toList().collect().ifEmpty("NONE").view()

        // temp_reads_ch.join(approved_samples_ch).view()

        // approved_samples_ch
        //     .join(reads_ch)
        //     .view()

    // emit:

}
