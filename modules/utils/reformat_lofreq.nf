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
process UTILS_REFORMAT_LOFREQ {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    container "ghcr.io/torch-consortium/magma/misc:2.0.0-alpha"

    input:
        tuple val(sampleName), path(lofreqVcf)

    output:
        tuple val(sampleName), path("*lofreq.reformat.corrected.vcf")

    script:
       
        """
        reformat_lofreq.py ${lofreqVcf} \\
            ${sampleName} \\
            ${sampleName}.lofreq.reformat.vcf

        reduce_strand_bias.py \\
            ${params.cutoff_strand_bias} \\
            ${sampleName}.lofreq.reformat.vcf  \\
            ${sampleName}.lofreq.reformat.corrected.vcf 
        """

    stub: 

        """
        touch ${sampleName}.lofreq.reformat.vcf
        touch ${sampleName}.lofreq.reformat.corrected.vcf 
        """ 

}
