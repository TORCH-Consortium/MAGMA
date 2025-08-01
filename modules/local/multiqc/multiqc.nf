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
process MULTIQC {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish
    label 'cpu_4_memory_16'

    input:
      path(multiqc_config)
        path("*")

    output:
        tuple path("multiqc_data"), path("multiqc_report.html")


    script:
        def config = multiqc_config ? "--config $multiqc_config" : ''

        def arg_skip_merge_analysis =  !params.skip_merge_analysis ? '' : '--skip_merge_analysis'

        """
        preprocess_multiqc_input.py   \\
            ${arg_skip_merge_analysis} \\
            --merged_cohort_stats joint.merged_cohort_stats.tsv \\
            --distance_matrix joint.ExDR.ExComplex.snp_dists.tsv

        ${params.multiqc_path} $config .
        """

    stub:
        """
        mkdir multiqc_data

        touch multiqc_report.html
        """
}
