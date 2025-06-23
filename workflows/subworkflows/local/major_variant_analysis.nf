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
include { TBPROFILER_VCF_PROFILE__COHORT } from "../../modules/local/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__COHORT)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__COHORT } from "../../modules/local/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__COHORT)

workflow MAJOR_VARIANT_ANALYSIS {
    take:
        merged_vcf_ch
        reformatted_lofreq_vcf_ch

    main:

        def resistanceDb =  []

        // merge_call_resistance
        TBPROFILER_VCF_PROFILE__COHORT(merged_vcf_ch, resistanceDb)
        TBPROFILER_COLLATE__COHORT(params.vcf_name, TBPROFILER_VCF_PROFILE__COHORT.out, resistanceDb)

        emit:
            major_variants_results_ch =  TBPROFILER_COLLATE__COHORT.out.per_sample_results
}
