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
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN7 } from "../../../modules/local/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN7 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN6 } from "../../../modules/local/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN6 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN5 } from "../../../modules/local/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN5 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN4 } from "../../../modules/local/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN4 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN3 } from "../../../modules/local/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN3 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN2 } from "../../../modules/local/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN2 )

include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN7 } from "../../../modules/local/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN7 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN6 } from "../../../modules/local/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN6 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN5 } from "../../../modules/local/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN5 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN4 } from "../../../modules/local/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN4 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN3 } from "../../../modules/local/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN3 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN2 } from "../../../modules/local/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN2 )

include { UTILS_SELECT_BEST_ANNOTATIONS } from "../../../modules/local/utils/select_best_annotations.nf" addParams ( params.UTILS_SELECT_BEST_ANNOTATIONS )

workflow OPTIMIZE_VARIANT_RECALIBRATION {
    take:
        analysisType
        select_variants_vcftuple_ch
        args_ch
        resources_files_ch
        resources_file_indexes_ch

    main:

        GATK_VARIANT_RECALIBRATOR__ANN7(analysisType,
                                    " -an DP -an AS_QD -an AS_MQ -an AS_FS -an AS_SOR -an AS_MQRankSum -an AS_ReadPosRankSum ",
                                    select_variants_vcftuple_ch,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )


        elim_ann7_ch = GATK_VARIANT_RECALIBRATOR__ANN7.out.annotationsLog
                        .join(GATK_VARIANT_RECALIBRATOR__ANN7.out.tranchesFile)


        UTILS_ELIMINATE_ANNOTATION__ANN7(analysisType, elim_ann7_ch)


//------------------------

        ann6_ch = UTILS_ELIMINATE_ANNOTATION__ANN7.out.reducedAnnotationsFile
                    .map { it.text }


        GATK_VARIANT_RECALIBRATOR__ANN6(analysisType,
                                    ann6_ch,
                                    select_variants_vcftuple_ch,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )


        elim_ann6_ch = GATK_VARIANT_RECALIBRATOR__ANN6.out.annotationsLog
                        .join(GATK_VARIANT_RECALIBRATOR__ANN6.out.tranchesFile)


        UTILS_ELIMINATE_ANNOTATION__ANN6(analysisType, elim_ann6_ch)


//------------------------

        ann5_ch = UTILS_ELIMINATE_ANNOTATION__ANN6.out.reducedAnnotationsFile
                    .map { it.text }


        GATK_VARIANT_RECALIBRATOR__ANN5(analysisType,
                                    ann5_ch,
                                    select_variants_vcftuple_ch,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )


        elim_ann5_ch = GATK_VARIANT_RECALIBRATOR__ANN5.out.annotationsLog
                        .join(GATK_VARIANT_RECALIBRATOR__ANN5.out.tranchesFile)


        UTILS_ELIMINATE_ANNOTATION__ANN5(analysisType, elim_ann5_ch)



//------------------------

        ann4_ch = UTILS_ELIMINATE_ANNOTATION__ANN5.out.reducedAnnotationsFile
                    .map { it.text }


        GATK_VARIANT_RECALIBRATOR__ANN4(analysisType,
                                    ann4_ch,
                                    select_variants_vcftuple_ch,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )

        elim_ann4_ch = GATK_VARIANT_RECALIBRATOR__ANN4.out.annotationsLog
                        .join(GATK_VARIANT_RECALIBRATOR__ANN4.out.tranchesFile)


        UTILS_ELIMINATE_ANNOTATION__ANN4(analysisType, elim_ann4_ch)

//------------------------

        ann3_ch = UTILS_ELIMINATE_ANNOTATION__ANN4.out.reducedAnnotationsFile
                    .map { it.text }


        GATK_VARIANT_RECALIBRATOR__ANN3(analysisType,
                                        ann3_ch,
                                        select_variants_vcftuple_ch,
                                        args_ch,
                                        resources_files_ch,
                                        resources_file_indexes_ch,
                                        params.ref_fasta,
                                        [params.ref_fasta_fai, params.ref_fasta_dict] )

        elim_ann3_ch = GATK_VARIANT_RECALIBRATOR__ANN3.out.annotationsLog
                        .join(GATK_VARIANT_RECALIBRATOR__ANN3.out.tranchesFile)


        UTILS_ELIMINATE_ANNOTATION__ANN3(analysisType, elim_ann3_ch)


//------------------------

        ann2_ch = UTILS_ELIMINATE_ANNOTATION__ANN3.out.reducedAnnotationsFile
                    .map { it.text }

        GATK_VARIANT_RECALIBRATOR__ANN2(analysisType,
                                    ann2_ch,
                                    select_variants_vcftuple_ch,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )

        elim_ann2_ch = GATK_VARIANT_RECALIBRATOR__ANN2.out.annotationsLog
                        .join(GATK_VARIANT_RECALIBRATOR__ANN2.out.tranchesFile)


        UTILS_ELIMINATE_ANNOTATION__ANN2(analysisType, elim_ann2_ch)


//------------------------
// NOTE: Choose the best set of annotations based on the highest minVQSLod score (closest to zero) for targetTruthSensitivity == 99.90
//------------------------

    recal_vcf_files_ch = GATK_VARIANT_RECALIBRATOR__ANN7.out.recalVcfTuple
                        .join(GATK_VARIANT_RECALIBRATOR__ANN6.out.recalVcfTuple)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN5.out.recalVcfTuple)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN4.out.recalVcfTuple)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN3.out.recalVcfTuple)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN2.out.recalVcfTuple)
                        .flatten()
                        .filter { it.class != String }
                        .collect()

    tranches_files_ch = GATK_VARIANT_RECALIBRATOR__ANN7.out.tranchesFile
                        .join(GATK_VARIANT_RECALIBRATOR__ANN6.out.tranchesFile)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN5.out.tranchesFile)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN4.out.tranchesFile)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN3.out.tranchesFile)
                        .join(GATK_VARIANT_RECALIBRATOR__ANN2.out.tranchesFile)
                        .flatten()
                        .filter { it.class != String }
                        .collect()


    annotations_and_tranches_json_files_ch = UTILS_ELIMINATE_ANNOTATION__ANN7.out.annotationsTranchesFile
                                .join(UTILS_ELIMINATE_ANNOTATION__ANN6.out.annotationsTranchesFile)
                                .join(UTILS_ELIMINATE_ANNOTATION__ANN5.out.annotationsTranchesFile)
                                .join(UTILS_ELIMINATE_ANNOTATION__ANN4.out.annotationsTranchesFile)
                                .join(UTILS_ELIMINATE_ANNOTATION__ANN3.out.annotationsTranchesFile)
                                .join(UTILS_ELIMINATE_ANNOTATION__ANN2.out.annotationsTranchesFile)
                                .flatten()
                                .filter { it.class != String }
                                .collect()




//------------------------

    UTILS_SELECT_BEST_ANNOTATIONS(params.vcf_name,
                                  annotations_and_tranches_json_files_ch,
                                  recal_vcf_files_ch,
                                  tranches_files_ch)


    //NOTE: We can run the other annotations process after the first ANN_7 process, in parallel but this is deffered to a future in interest of Engg. effort.
    emit:
        optimized_vqsr_ch = select_variants_vcftuple_ch
                            .join(UTILS_SELECT_BEST_ANNOTATIONS.out.bestRecalVcfTuple)
                            .join(UTILS_SELECT_BEST_ANNOTATIONS.out.bestTranchesFile)


}
