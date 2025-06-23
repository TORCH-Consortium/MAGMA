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
include { GATK_SELECT_VARIANTS as GATK_SELECT_VARIANTS__PHYLOGENY } from "../../modules/local/gatk/select_variants.nf" addParams( params.GATK_SELECT_VARIANTS__PHYLOGENY )
include { GATK_VARIANTS_TO_TABLE } from "../../modules/local/gatk/variants_to_table.nf" addParams( params.GATK_VARIANTS_TO_TABLE )
include { SNPSITES } from "../../modules/local/snpsites/snpsites.nf" addParams( params.SNPSITES )
include { SNPDISTS } from "../../modules/local/snpdists/snpdists.nf" addParams( params.SNPDISTS )
include { IQTREE } from "../../modules/local/iqtree/iqtree.nf" addParams( params.IQTREE )

workflow PHYLOGENY_ANALYSIS {

    take:
        prefix_ch
        arg_files_ch
        vcf_ch

    main:

        args_ch = arg_files_ch
            .filter {  (it.getExtension()  == "gz") || (it.getExtension()  == "list") }
            .reduce (" ") { a, b -> {
                    if (b.class == sun.nio.fs.UnixPath) {
                        return "$a -XL ${b.getName()} "

                    } else {
                        return ""
                    }
                }
            }
            //.dump(tag: "PHYLOGENY_ANALYSIS args_ch: ", pretty: true)


        resources_files_ch = arg_files_ch
            .filter {  (it.getExtension()  == "gz") || (it.getExtension()  == "list") }
            .collect()
            .ifEmpty([])
            //.dump(tag: "PHYLOGENY_ANALYSIS resources_files_ch: ", pretty: true)


        resources_file_indexes_ch = arg_files_ch
            .filter {  it.getExtension()  == "tbi" }
            .collect()
            .ifEmpty([])
            //.dump(tag: "PHYLOGENY_ANALYSIS resources_file_indexes_ch: ", pretty: true)


        // merge_phylogeny_prep_inccomplex
        GATK_SELECT_VARIANTS__PHYLOGENY('SNP',
                            prefix_ch,
                            vcf_ch,
                            args_ch,
                            resources_files_ch,
                            resources_file_indexes_ch,
                            params.ref_fasta,
                            [params.ref_fasta_fai, params.ref_fasta_dict])

        GATK_VARIANTS_TO_TABLE(prefix_ch, GATK_SELECT_VARIANTS__PHYLOGENY.out)

        SNPSITES(prefix_ch, GATK_VARIANTS_TO_TABLE.out)

        SNPDISTS(prefix_ch, SNPSITES.out)

        // merge_iqtree_inccomplex
        IQTREE(prefix_ch, SNPSITES.out)


    emit:
        snpsites_tree_tuple = SNPSITES.out.join(IQTREE.out.tree_tuple)
  snp_dists_ch = SNPDISTS.out.snp_dists_file
}
