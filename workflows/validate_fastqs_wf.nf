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
include { FASTQ_VALIDATOR } from '../modules/local/fastq_utils/validator.nf' addParams ( params.FASTQ_VALIDATOR  )
include { UTILS_FASTQ_COHORT_VALIDATION } from '../modules/local/utils/fastq_cohort_validation.nf' addParams ( params.UTILS_FASTQ_COHORT_VALIDATION  )

workflow VALIDATE_FASTQS_WF {
    take:
         samplesheet_json
         ready

    main:

    //NOTE: Expected structure of input CSV samplesheet, including all optional fields
    //   0     1       2       3    4  5     6      7       8
    // Study,Sample,Library,Attempt,R1,R2,Flowcell,Lane,Index Sequence



    fastqs_ch = samplesheet_json
        .splitJson()
        .map {
            if (it.R2) {
                [it.magma_sample_name, [it.R1, it.R2]]
            } else {

                [it.magma_sample_name, [it.R1]]
            }
        }.transpose()




    FASTQ_VALIDATOR( fastqs_ch, ready )


    UTILS_FASTQ_COHORT_VALIDATION( FASTQ_VALIDATOR.out.fastq_report.collect(), samplesheet_json )


    validated_fastqs_ch = UTILS_FASTQ_COHORT_VALIDATION.out.magma_analysis_json
                            .splitJson()
                            .filter {it.value.fastqs_approved}
                            .map {
                                    if (it.value.R2) {
                [it.value.magma_sample_name, [bam_rg_string: it.value.magma_bam_rg_string, paired: true], [it.value.R1, it.value.R2]]
                                    } else {
                [it.value.magma_sample_name, [bam_rg_string: it.value.magma_bam_rg_string, paired: false] , [it.value.R1]]
                                    }
                            }

     emit:

        approved_fastqs_ch = validated_fastqs_ch


}
