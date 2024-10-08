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
include { BWA_MEM } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM)

workflow MAP_WF {
    take:
        samples_ch

    main:

        BWA_MEM(samples_ch,
                  params.ref_fasta,
                  [params.ref_fasta_dict,
                   params.ref_fasta_amb,
                   params.ref_fasta_ann,
                   params.ref_fasta_bwt,
                   params.ref_fasta_fai,
                   params.ref_fasta_pac,
                   params.ref_fasta_sa])


    emit:
        sorted_reads_ch = BWA_MEM.out

}
