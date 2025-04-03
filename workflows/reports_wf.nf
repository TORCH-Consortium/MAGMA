/*
 * Copyright (c) 2021-2024 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
 *
 * This file is part of MAGMA pipeline, see https://github.com/TORCH-Consortium/MAGMA
 *
 * For quick overview of GPL-3 license, please refer
 * https://www.tldrlegal.com/license/gnu-general-public-license-v3-gpl-3
 *
 * - You MUST keep this license with original authors in your copy
 * - You MUST acknowledge the original source of this software
 * - You MUST state significant changes made to the original software
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program . If not, see <http://www.gnu.org/licenses/>.
 */
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
