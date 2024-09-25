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
include { GATK_SELECT_VARIANTS as  GATK_SELECT_VARIANTS__SNP } from "../../modules/gatk/select_variants.nf" addParams ( params.GATK_SELECT_VARIANTS__SNP )
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR__SNP } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__SNP )
include { GATK_APPLY_VQSR as GATK_APPLY_VQSR__SNP } from "../../modules/gatk/apply_vqsr.nf" addParams ( params.GATK_APPLY_VQSR__SNP )
include { GATK_SELECT_VARIANTS__EXCLUSION as  GATK_SELECT_VARIANTS__EXCLUSION__SNP } from "../../modules/gatk/select_variants__exclusion.nf" addParams ( params.GATK_SELECT_VARIANTS__EXCLUSION__SNP )
include { OPTIMIZE_VARIANT_RECALIBRATION } from "./optimize_variant_recalibration.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__SNP )

workflow SNP_ANALYSIS {

    take:
        cohort_vcf_and_index_ch

    main:

        // merge_select_snp
        GATK_SELECT_VARIANTS__SNP('SNP',
                                'raw',
                                cohort_vcf_and_index_ch,
                                "",
                                [],
                                [],
                                params.ref_fasta,
                                [params.ref_fasta_fai, params.ref_fasta_dict] )


        // merge_vqsr_snp

        arg_files_ch = Channel.of(["coll2014,known=false,training=true,truth=true,prior=15.0", file(params.coll2014_vcf), file(params.coll2014_vcf_tbi)],
                              ["coll2018,known=false,training=true,truth=true,prior=15.0", file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)],
                              ["Napier2020,known=false,training=true,truth=true,prior=15.0", file(params.napier2020_vcf), file(params.napier2020_vcf_tbi)],
                              ["Benavente2015,known=true,training=false,truth=false,prior=5.0", file(params.benavente2015_vcf), file(params.benavente2015_vcf_tbi)])
            .ifEmpty([])
            .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1], it[2] ] : [] }
            .flatten()
            .dump(tag:"SNP_ANALYSIS arg_files_ch : ", pretty:true)


        args_ch = arg_files_ch
            .filter { it.class == org.codehaus.groovy.runtime.GStringImpl }
            .reduce { a, b -> "$a --resource:$b " }
            .ifEmpty("")
            .dump(tag:"SNP_ANALYSIS args_ch : ", pretty:true)


        resources_files_ch = arg_files_ch
            .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
            .filter {  it.getExtension()  == "gz" }
            .collect()
            .ifEmpty([])
            .dump(tag:"SNP_ANALYSIS resources_files_ch : ", pretty:true)

        resources_file_indexes_ch = arg_files_ch
            .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
            .filter {  it.getExtension()  == "tbi" }
            .collect()
            .ifEmpty([])
            .dump(tag:"SNP_ANALYSIS resources_file_indexes_ch : ", pretty:true)



        if(!params.skip_variant_recalibration) {

            OPTIMIZE_VARIANT_RECALIBRATION('SNP',
                                        GATK_SELECT_VARIANTS__SNP.out.variantsVcfTuple,
                                        args_ch,
                                        resources_files_ch,
                                        resources_file_indexes_ch)

            vqsr_ch = OPTIMIZE_VARIANT_RECALIBRATION.out.optimized_vqsr_ch

        } else {

            GATK_VARIANT_RECALIBRATOR__SNP('SNP',
                                    " -an DP -an AS_MQ ",
                                        GATK_SELECT_VARIANTS__SNP.out.variantsVcfTuple,
                                        args_ch,
                                        resources_files_ch,
                                        resources_file_indexes_ch,
                                        params.ref_fasta,
                                        [params.ref_fasta_fai, params.ref_fasta_dict] )

            vqsr_ch = GATK_SELECT_VARIANTS__SNP.out.variantsVcfTuple
                .join(GATK_VARIANT_RECALIBRATOR__SNP.out.recalVcfTuple)
                .join(GATK_VARIANT_RECALIBRATOR__SNP.out.tranchesFile)

        }


        // merge_apply_vqsr_snp
        GATK_APPLY_VQSR__SNP('SNP',
                            vqsr_ch,
                            params.ref_fasta,
                            [params.ref_fasta_fai, params.ref_fasta_dict])

        GATK_SELECT_VARIANTS__EXCLUSION__SNP('SNP',
                                            GATK_APPLY_VQSR__SNP.out.filteredVcfTuple,
                                            params.rrna_list,
                                            params.ref_fasta,
                                            [params.ref_fasta_fai, params.ref_fasta_dict])


    emit:
        snp_exc_vcf_ch = GATK_SELECT_VARIANTS__EXCLUSION__SNP.out
        snp_inc_vcf_ch = GATK_APPLY_VQSR__SNP.out.filteredVcfTuple

}
