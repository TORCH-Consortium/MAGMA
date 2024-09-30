nextflow.enable.dsl = 2

include { MULTIQC } from '../modules/multiqc/multiqc.nf' addParams (params.MULTIQC)
include { UTILS_SUMMARIZE_RESISTANCE_RESULTS } from '../modules/utils/summarize_resistance_results.nf' addParams (params.UTILS_SUMMARIZE_RESISTANCE_RESULTS)
include { UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION } from "../modules/utils/summarize_resistance_results_mixed_infection.nf" addParams (params.UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION)


workflow REPORTS_WF {
    take:
         reports_fastqc_ch
         merged_cohort_stats_ch
         major_variants_results_ch
         minor_variants_results_ch
         structural_variants_results_ch
	 snp_distances_ch

    main:
	
	ch_multiqc_files = Channel.empty()

	UTILS_SUMMARIZE_RESISTANCE_RESULTS(
            merged_cohort_stats_ch,
            major_variants_results_ch,
            minor_variants_results_ch,
            structural_variants_results_ch
        )

        UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION(
            merged_cohort_stats_ch,
            minor_variants_results_ch,
            structural_variants_results_ch
        )
        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

        ch_multiqc_files = ch_multiqc_files.mix(
	reports_fastqc_ch,
	merged_cohort_stats_ch,
	snp_distances_ch
	)

        MULTIQC(ch_multiqc_config,
	ch_multiqc_files.collect())
}
