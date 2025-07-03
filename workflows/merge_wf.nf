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
include { PREPARE_COHORT_VCF } from "../subworkflows/local/prepare_cohort_vcf.nf"
include { SNP_ANALYSIS } from "../subworkflows/local/snp_analysis.nf"
include { INDEL_ANALYSIS } from "../subworkflows/local/indel_analysis.nf"
include { GATK_MERGE_VCFS as GATK_MERGE_VCFS__INC } from "../modules/local/gatk/merge_vcfs.nf" addParams ( params.GATK_MERGE_VCFS )
include { MAJOR_VARIANT_ANALYSIS } from "../subworkflows/local/major_variant_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS__INCCOMPLEX } from "../subworkflows/local/phylogeny_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS__EXCOMPLEX } from "../subworkflows/local/phylogeny_analysis.nf"
include { CLUSTER_ANALYSIS as CLUSTER_ANALYSIS__INCCOMPLEX } from "../subworkflows/local/cluster_analysis.nf"
include { CLUSTER_ANALYSIS as  CLUSTER_ANALYSIS__EXCOMPLEX } from "../subworkflows/local/cluster_analysis.nf"


workflow MERGE_WF {
    take:
        gvcf_ch
        reformatted_lofreq_vcfs_tuple_ch
        approved_samples_ch

    main:

        //---------------------------------------------------------------------------------
        // Filter the approved samples
        //---------------------------------------------------------------------------------

        //NOTE: Reshape the flattened output of gvch_ch into the tuples of [sampleName, gvcf, gvcf.tbi]
        collated_gvcfs_ch = gvcf_ch
                                .flatten()
                                .collate(3)
                                .dump(tag:'MERGE_WF: collated_gvcfs_ch', pretty: true)
                                //.collectFile(name: "$params.outdir/collated_gvcfs_ch.txt")

        //NOTE: Join the fully approved samples with the gvcf channel
        selected_gvcfs_ch = collated_gvcfs_ch.join(approved_samples_ch)
                                        .flatten()
                                        .dump(tag:'MERGE_WF: selected_gvcfs_ch', pretty: true)


        //NOTE: Filter only file type values and send to MERGE_WF
        //FIXME refactor the filtering logic NOT to rely upon the exact classnames
        filtered_selected_gvcfs_ch = selected_gvcfs_ch
                                        .filter { it.class != String }
                                        .collect()
                                        //.dump(tag:'MERGE_WF: filtered_selected_gvcfs_ch', pretty: true)
                                        //.collectFile(name: "$params.outdir/selected_gvcfs_ch")

        //filtered_selected_gvcfs_ch.view()


       //---------------------------------------------------------------------------------

        PREPARE_COHORT_VCF( filtered_selected_gvcfs_ch )

        SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        INDEL_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

    if(!params.skip_variant_recalibration )  {

        //Use the recalibrated SNPs
        merge_vcf_ch = SNP_ANALYSIS.out.snp_inc_vcf_ch
                            .join(INDEL_ANALYSIS.out.indel_vcf_ch)
                            .dump(tag:'MERGE_WF: merge_inc_vcf_ch : ', pretty: true)

        snp_exc_vcf_ch = SNP_ANALYSIS.out.snp_exc_vcf_ch


    } else {

        //Do NOT use the recalibrated SNPs
        merge_vcf_ch = SNP_ANALYSIS.out.snp_vcf_ch
                            .join(INDEL_ANALYSIS.out.indel_vcf_ch)
                            .dump(tag:'MERGE_WF: merge_inc_vcf_ch : ', pretty: true)

        snp_exc_vcf_ch = SNP_ANALYSIS.out.snp_vcf_ch

    }
        // merge_snp_indel_vcf
        GATK_MERGE_VCFS__INC(merge_vcf_ch)

        MAJOR_VARIANT_ANALYSIS(GATK_MERGE_VCFS__INC.out, reformatted_lofreq_vcfs_tuple_ch)


        //----------
        // Exclude complex regions for downstream analysis
        //----------


            if (!params.skip_phylogeny_and_clustering) {

                excomplex_prefix_ch = Channel.of('ExDR.ExComplex')

                excomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)],
                                                            [file(params.excluded_loci_list)])
                                                        .ifEmpty([])
                                                        .flatten()
                                                        //.dump(tag:'MERGE_WF: excomplex_exclude_interval_ref_ch', pretty: true)


                // excomplex_exclude_interval_ref_ch.view{ it -> "\n\n MAGMA-LOG MERGE_WF excomplex_exclude_interval_ref_ch: $it \n\n"}




                PHYLOGENY_ANALYSIS__EXCOMPLEX(excomplex_prefix_ch,
                                              excomplex_exclude_interval_ref_ch,
                                              snp_exc_vcf_ch)



                CLUSTER_ANALYSIS__EXCOMPLEX(PHYLOGENY_ANALYSIS__EXCOMPLEX.out.snpsites_tree_tuple, excomplex_prefix_ch)

        }



    //------------------------------------------------------
    //     Include complex regions for downstream analysis
    // NOTE: This is officially recommended as per XBS paper
    //------------------------------------------------------

        if(!params.skip_complex_regions) {

            if (!params.skip_phylogeny_and_clustering) {

                inccomplex_prefix_ch = Channel.of('ExDR.IncComplex')

                inccomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)])
                                                        .ifEmpty([])
                                                        .flatten()




                //NOTE: Both phylogenies should be excluding DR and excluding rRNA, then it is again filtered in two datasets one including complex regions and one excluding complex regions.
                //Ergo PHYLOGENY_...__INCCOMPLEX should take snp_exc_vcf_ch. Refer https://github.com/TORCH-Consortium/MAGMA/pull/114#discussion_r947732253
                PHYLOGENY_ANALYSIS__INCCOMPLEX(inccomplex_prefix_ch,
                                               inccomplex_exclude_interval_ref_ch,
                                               snp_exc_vcf_ch)

                CLUSTER_ANALYSIS__INCCOMPLEX(PHYLOGENY_ANALYSIS__INCCOMPLEX.out.snpsites_tree_tuple, inccomplex_prefix_ch)
            }

        }



    emit:
        major_variants_results_ch =  MAJOR_VARIANT_ANALYSIS.out.major_variants_results_ch
        snps_dists_ch = PHYLOGENY_ANALYSIS__EXCOMPLEX.out.snp_dists_ch
}
