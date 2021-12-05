include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR_ANN7;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR_ANN6;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR_ANN5;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR_ANN4;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR_ANN3;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR_ANN2
} from "../../modules/gatk/variant_recalibrator.nf"


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



    emit:

        optimized_vqsr_ch = select_variants_vcftuple_ch
                            .join(GATK_VARIANT_RECALIBRATOR__ANN7.out.recalVcfTuple)
                            .join(GATK_VARIANT_RECALIBRATOR__ANN7.out.tranchesFile)


}
