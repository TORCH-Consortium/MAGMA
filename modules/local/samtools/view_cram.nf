/*
 * Copyright (c) 2021-2026 MAGMA pipeline authors, see https://doi.org/10.1371/journal.pcbi.1011648
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

// SAMTOOLS_VIEW_CRAM — optional BAM → CRAM checkpoint (Pattern 1 of the
// BAM-checkpoints spec, abc-universe/specs/active/magma-bam-checkpoints.md).
//
// Drop-in replacement for SAMTOOLS_INDEX at the post-MarkDup-or-BQSR
// junction in call_wf.nf. Emits the identical channel shape
// `tuple(sampleName, index, alignment)` so every downstream consumer
// (HaplotypeCaller, minor-variants HC) is agnostic to whether the input
// is BAM+BAI or CRAM+CRAI. GATK reads CRAM directly via `-I <cram> -R <ref>`.
//
// Lossless by default (cram_version=3.1, default codec set). Quality-score
// binning is opt-in via params.cram_lossy_qualities (use only for archival
// of already-called cohorts; not safe for re-calling or re-BQSR).

process SAMTOOLS_VIEW_CRAM {
    tag "${sampleName}"
    label 'cpu_4_memory_4'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(bam)
        path(ref_fasta)
        path("*")                    // ref_fasta_fai etc. (mirrors HaplotypeCaller's idiom)

    output:
        tuple val(sampleName), path("*.cram.crai"), path("*.cram")

    script:
        // Lossless CRAM 3.1 by default. With `cram_lossy_qualities=true`,
        // bin quality scores to 8 bins (samtools' --output-fmt-option
        // lossy_names=0,store_md=1,store_nm=1,...) — smaller but lossy.
        def lossy_opts = params.cram_lossy_qualities
            ? '--output-fmt-option lossy_names=0,store_md=1,store_nm=1,seqs_per_slice=10000'
            : ''
        """
        ${params.samtools_path} view \\
            -C \\
            -T ${ref_fasta} \\
            --threads ${task.cpus} \\
            --output-fmt cram,version=3.1 \\
            ${lossy_opts} \\
            -o ${sampleName}.cram \\
            ${bam}

        ${params.samtools_path} index \\
            -@ ${task.cpus} \\
            ${sampleName}.cram
        """

    stub:
        """
        echo "samtools view -C -T ${ref_fasta} -o ${sampleName}.cram ${bam}"
        echo "samtools index ${sampleName}.cram"
        touch ${sampleName}.cram
        touch ${sampleName}.cram.crai
        """
}
