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
process TBPROFILER_FASTQ_PROFILE {
    tag "$sampleName"
    label 'process_medium'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish, pattern: "!*json"

    input:
        tuple val(sampleName), val(meta), path(reads)

    output:
        tuple val(meta), path("results/*.json"), emit: json
        tuple val(meta), path("results/*.csv") , emit: csv, optional: true
        tuple val(meta), path("results/*.txt") , emit: txt, optional: true

        //NOTE: Disable these outputs since we are not using them
        // path "versions.yml"                    , emit: versions
        // tuple val(meta), path("bam/*.bam")     , emit: bam
        // tuple val(meta), path("vcf/*.vcf.gz")  , emit: vcf

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args   ?: ''
        prefix   = task.ext.prefix ?: "${sampleName}"
        def input_reads = meta.single_end ? "--read1 $reads" : "--read1 ${reads[0]} --read2 ${reads[1]}"
        """
        tb-profiler \\
            profile \\
            $args \\
            ${params.arguments} \\
            --prefix ${prefix} \\
            --threads $task.cpus \\
            $input_reads

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            tbprofiler:  \$( echo \$(tb-profiler --version 2>&1) | sed 's/TBProfiler version //')
        END_VERSIONS
        """

    stub:
        """
        mkdir results
        touch results/${sampleName}.results.json
        """

}
