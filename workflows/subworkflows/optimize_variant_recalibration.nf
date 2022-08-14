include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN7 } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN7 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN6 } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN6 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN5 } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN5 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN4 } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN4 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN3 } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN3 )
include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN2 } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__ANN2 )

include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN7 } from "../../modules/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN7 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN6 } from "../../modules/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN6 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN5 } from "../../modules/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN5 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN4 } from "../../modules/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN4 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN3 } from "../../modules/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN3 )
include { UTILS_ELIMINATE_ANNOTATION as  UTILS_ELIMINATE_ANNOTATION__ANN2 } from "../../modules/utils/eliminate_annotation.nf" addParams ( params.UTILS_ELIMINATE_ANNOTATION__ANN2 )


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


//------------------------

        UTILS_ELIMINATE_ANNOTATION__ANN7(params.vcf_name,
                                         analysisType,
                                         GATK_VARIANT_RECALIBRATOR__ANN7.out.annotationsLog,
                                         GATK_VARIANT_RECALIBRATOR__ANN7.out.tranchesFile )


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

//------------------------

        UTILS_ELIMINATE_ANNOTATION__ANN6(params.vcf_name,
                                         analysisType,
                                         GATK_VARIANT_RECALIBRATOR__ANN6.out.annotationsLog,
                                         GATK_VARIANT_RECALIBRATOR__ANN6.out.tranchesFile)


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



//------------------------

        UTILS_ELIMINATE_ANNOTATION__ANN5(params.vcf_name,
                                        analysisType,
                                        GATK_VARIANT_RECALIBRATOR__ANN5.out.annotationsLog,
                                        GATK_VARIANT_RECALIBRATOR__ANN5.out.tranchesFile)


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

//------------------------

        UTILS_ELIMINATE_ANNOTATION__ANN4(params.vcf_name,
                                        analysisType,
                                        GATK_VARIANT_RECALIBRATOR__ANN4.out.annotationsLog,
                                        GATK_VARIANT_RECALIBRATOR__ANN4.out.tranchesFile)

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

//------------------------


        UTILS_ELIMINATE_ANNOTATION__ANN3(params.vcf_name,
                                        analysisType,
                                        GATK_VARIANT_RECALIBRATOR__ANN3.out.annotationsLog,
                                        GATK_VARIANT_RECALIBRATOR__ANN3.out.tranchesFile)

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

//------------------------

        UTILS_ELIMINATE_ANNOTATION__ANN2(params.vcf_name,
                                        analysisType,
                                        GATK_VARIANT_RECALIBRATOR__ANN2.out.annotationsLog,
                                        GATK_VARIANT_RECALIBRATOR__ANN2.out.tranchesFile)



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
