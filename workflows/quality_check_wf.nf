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
include { FASTQC              } from '../modules/fastqc/fastqc.nf' addParams (params.FASTQC)
include { NTMPROFILER_PROFILE } from '../modules/ntmprofiler/profile.nf' addParams (params.NTMPROFILER_PROFILE)
include { NTMPROFILER_COLLATE } from '../modules/ntmprofiler/collate.nf' addParams (params.NTMPROFILER_COLLATE)

include { TBPROFILER_FASTQ_PROFILE as EXP_TBPROFILER_FASTQ_PROFILE } from '../modules/tbprofiler/fastq_profile.nf' addParams (params.EXP_TBPROFILER_FASTQ_PROFILE)
include { TBPROFILER_COLLATE as EXP_TBPROFILER_FASTQ_COLLATE } from '../modules/tbprofiler/collate.nf' addParams (params.EXP_TBPROFILER_FASTQ_COLLATE)

include { SPOTYPING as EXP_SPOTYPING } from '../modules/spotyping/main.nf' addParams (params.EXP_SPOTYPING)
include { RDANALYZER as EXP_RDANALYZER } from '../modules/rdanalyzer/main.nf' addParams (params.EXP_RDANALYZER)

workflow QUALITY_CHECK_WF {

    take:
        reads_ch

    main:

        FASTQC(reads_ch)

        NTMPROFILER_PROFILE( reads_ch )

        NTMPROFILER_COLLATE( params.vcf_name,
                             NTMPROFILER_PROFILE.out.profile_json.collect() )

        if (!params.skip_tbprofiler_fastq) {

            EXP_TBPROFILER_FASTQ_PROFILE( reads_ch )

            EXP_TBPROFILER_FASTQ_COLLATE( params.vcf_name,
                                          EXP_TBPROFILER_FASTQ_PROFILE.out.profile_json.collect(),
                                          [] )
        }


        if(!params.exp_skip_rdanalyzer) {
            EXP_RDANALYZER( reads_ch )
        }


        if(!params.exp_skip_spotyping) {
            EXP_SPOTYPING( reads_ch )
        }


    emit:
        reports_fastqc_ch =  FASTQC.out.collect()

}
