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
process BWA_MEM {
    tag "${sampleName}"
    label 'cpu_8_memory_8'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), val(meta), path(sampleReads)
        path(reference)
        path("*")

    output:
        tuple val(sampleName), path("*.sorted_reads.bam")


    script:

        """
        ${params.bwa_path} mem \\
            -M \\
            ${params.arguments} \\
            -t ${task.cpus} \\
            -R "${meta.bam_rg_string}"  \\
            ${reference} \\
            ${sampleReads} \\
        | ${params.samtools_path} sort \\
            -@ ${task.cpus} \\
            -O BAM \\
            -o ${sampleName}.sorted_reads.bam -
        """

    stub:

        """
        echo "${params.bwa_path} mem \\
            -M \\
            -t ${task.cpus} \\
            -R ${meta.bam_rg_string}  \\
            ${reference} \\
            ${sampleReads} \\
        | ${params.samtools_path} sort \\
            -@ ${task.cpus} \\
            -O BAM \\
            -o ${sampleName}.sorted_reads.bam -
        "

        touch ${sampleName}.sorted_reads.bam
        """

}
