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
include { GATK_COMBINE_GVCFS } from "../../modules/local/gatk/combine_gvcfs.nf" addParams ( params.GATK_COMBINE_GVCFS )
include { GATK_GENOTYPE_GVCFS } from "../../modules/local/gatk/genotype_gvcfs.nf" addParams ( params.GATK_GENOTYPE_GVCFS )
include { SNPEFF } from "../../modules/local/snpeff/snpeff.nf" addParams ( params.SNPEFF )
include { BGZIP } from "../../modules/local/bgzip/bgzip.nf" addParams( params.BGZIP )
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE__COHORT } from "../../modules/local/gatk/index_feature_file" addParams( params.GATK_INDEX_FEATURE_FILE__COHORT )


workflow PREPARE_COHORT_VCF {
    take:
        cohort_gvcfs_ch


    main:
        //FIXME save the string to an intermediate file
        gvcfs_string_ch = cohort_gvcfs_ch
            .flatten()
            .filter {  file(it).getExtension()  == "gz" }
            .map { it -> file(it).name }
            .reduce { a, b -> "$a --variant $b " }
            .dump(tag: "PREPARE_COHORT_VCF gvcfs_string_ch:", pretty: true)


        if (params.use_ref_gvcf) {
            refExitRifGvcf = file(params.ref_gvcf, checkIfExists: true)
        } else {
            refExitRifGvcf = []
        }

        if (params.use_ref_gvcf) {
            refExitRifGvcfTbi = file(params.ref_gvcf_tbi, checkIfExists: true)
        } else {
            refExitRifGvcfTbi = []
        }

        // merge_combine
        GATK_COMBINE_GVCFS(params.vcf_name,
                        gvcfs_string_ch,
                        cohort_gvcfs_ch,
                        params.ref_fasta,
                        refExitRifGvcf,
                        refExitRifGvcfTbi,
                        [params.ref_fasta_fai, params.ref_fasta_dict])


        // merge_genotype
        GATK_GENOTYPE_GVCFS(GATK_COMBINE_GVCFS.out, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])

        // merge_snpeff_annotate
        SNPEFF(GATK_GENOTYPE_GVCFS.out, params.ref_fasta)
        BGZIP(SNPEFF.out)
        GATK_INDEX_FEATURE_FILE__COHORT(BGZIP.out, '')

    emit:
        cohort_vcf_and_index_ch = GATK_INDEX_FEATURE_FILE__COHORT.out.sample_vcf_tuple
}
