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
process CLUSTERPICKER {
    tag "snpCount - ${snpCount}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(joint_name), path(fasta), path(newickTree)
        val(snpCount)
        val(prefix)

    output:
        tuple val(joint_name), path("*${snpCount}SNPcluster*")

    shell:
        //NOTE: The java_opts are propogated via bash as per the wrapper script in bioconda
        // https://github.com/bioconda/bioconda-recipes/blob/master/recipes/clusterpicker/cluster-picker.sh


        //NOTE: As per the instructions in the cluster-picker manual
        // cluster-picker input-fasta.fas                                    input-tree.nwk                                  bootstrap bootstrap genetic-distance max_cluster_size dist
        // cluster-picker $JOINT_NAME.95X.variable.IncComplex.5SNPcluster.fa $JOINT_NAME.95X.IncComplex.5SNPcluster.treefile 0         0         $FRACTION1       0              gap

        '''

        cp !{fasta} !{fasta.getBaseName()}.!{snpCount}SNPcluster.fa

        cp !{newickTree} !{newickTree.getBaseName()}.!{snpCount}SNPcluster.treefile

        # computes the genetic-distance argument
        FRACTION=$(echo !{snpCount}/$(cat !{fasta} | head -n 2 | tail -n 1| wc -c) | bc -l)

        jvm_mem_opts="-Xmx!{task.memory.giga}G"

        !{params.clusterpicker_path} \\
            !{fasta} \\
            !{newickTree} \\
            !{params.bootstrap_1} \\
            !{params.bootstrap_2} \\
            $FRACTION \\
            !{params.max_cluster_size} \\
            !{params.algorithm}

        cp !{joint_name}.!{prefix}_clusterPicks_log.txt !{joint_name}.!{prefix}_!{snpCount}SNPclusterPicks_log.txt
        cp !{joint_name}.!{prefix}_clusterPicks.nwk !{joint_name}.!{prefix}_!{snpCount}SNPclusterPicks.nwk
        cp !{joint_name}.!{prefix}_clusterPicks.nwk.figTree !{joint_name}.!{prefix}_!{snpCount}SNPclusterPicks.nwk.figTree
        cp !{joint_name}.variable.!{prefix}.fa_!{joint_name}.!{prefix}_clusterPicks.fas !{joint_name}.variable.!{prefix}.fa_!{joint_name}.!{prefix}_!{snpCount}SNPclusterPicks.fas

        '''

    stub:

        """
        """

}
