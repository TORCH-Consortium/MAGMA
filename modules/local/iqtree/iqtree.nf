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
process IQTREE {
    label 'cpu_8_memory_16'
    tag "${joint_name}"

    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:

        val(prefix)
        tuple val(joint_name), path(fasta)

    output:

        tuple val(joint_name), path("${joint_name}.${prefix}.treefile"), emit: tree_tuple
        path("${joint_name}.${prefix}.bionj")
        path("${joint_name}.${prefix}.ckp.gz")
        path("${joint_name}.${prefix}.iqtree")
        path("${joint_name}.${prefix}.log")
        path("${joint_name}.${prefix}.mldist")
        path("${joint_name}.${prefix}.model.gz")
        path("${joint_name}.${prefix}.treefile")


    script:
        if(params.iqtree_custom_argument) {
            arguments = params.iqtree_custom_argument
        } else if(params.iqtree_standard_bootstrap) {
            arguments = '-b 1000'
        } else if(params.iqtree_fast_ml_only) {
            arguments = '-fast'
        } else if(params.iqtree_fast_bootstrapped_phylogeny) {
            arguments = '-bb 1000 -alrt 1000'
        } else if(params.iqtree_accurate_ml_only) {
            arguments = '-allnni'
        } else {
        //NOTE: Use iqtree_accurate_ml_only as the default
            arguments = '-allnni'
        }

        """
        ${params.iqtree_path} \\
            -s ${fasta} \\
            -T AUTO \\
            ${arguments} \\
            --prefix ${joint_name}.${prefix}
        """

    stub:

        """
        touch ${joint_name}.${prefix}.bionj
        touch ${joint_name}.${prefix}.ckp.gz
        touch ${joint_name}.${prefix}.iqtree
        touch ${joint_name}.${prefix}.log
        touch ${joint_name}.${prefix}.mldist
        touch ${joint_name}.${prefix}.model.gz
        touch ${joint_name}.${prefix}.treefil
        """

}
