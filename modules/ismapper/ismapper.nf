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
process ISMAPPER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish


    input:
        tuple val(sampleName), val(meta), path(reads)
        path(ref_gbk)
        path(ref_fasta)
        path(query_multifasta)


    output:
        tuple val(sampleName), path("*.ismapper.vcf"), emit: formatted_vcf
        tuple val(sampleName), path("ismapper"), emit: results_dir

    script:

        """

        ismap \\
            --t $task.cpus \\
            --output_dir $sampleName \\
            --queries $query_multifasta \\
            --reference $ref_gbk \\
            --reads $reads

        # Reorganize output files
        mkdir ismapper
        mv $sampleName/*/* ismapper/


        # Reformat the is-mapper output to a standard VCF
        # FIXME the IS element name is hard-coded
        convert_ismapper_to_vcf.py \\
            --is_mapper_dir ismapper/IS6110/ \\
            --reference_file $ref_fasta \\
            --query_file $query_multifasta \\
            --output_vcf_file ${sampleName}.ismapper.vcf

        """
}
