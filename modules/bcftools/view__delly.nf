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
process BCFTOOLS_VIEW__DELLY {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), path(vcf)

    output:
        tuple val(sampleName), path("*.vcf.gz"), path("*.vcf.gz.csi")

    shell:

        '''
        !{params.bcftools_path} view !{params.arguments}  !{sampleName}.delly.vcf -o !{sampleName}.filtered.delly.vcf
        bgzip !{sampleName}.filtered.delly.vcf
        !{params.bcftools_path} index !{sampleName}.filtered.delly.vcf.gz
        '''

    stub:

        """
        touch ${sampleName}.vcf.gz
        touch ${sampleName}.vcf.gz.csi
        """

}
