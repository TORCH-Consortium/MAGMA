include { GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN7;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN6;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN5;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN4;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN3;
          GATK_VARIANT_RECALIBRATOR as  GATK_VARIANT_RECALIBRATOR__ANN2
} from "../../modules/gatk/variant_recalibrator.nf"


    //TODO: Perhaps better to place this in the util functions
        def eliminateLeastInformativeAnnotation(logFile) {
                logFile.eachLine {
                        if (it.contains("VariantDataManager - Annotation order is")) {
                            def chunkedString = it.split(":")
                            def orderedAnnotationsArray = chunkedString[chunkedString.size() - 1]
                            def orderedAnnotationsString = orderedAnnotationsArray.replace('[', ' ').replace(']', ' ')

                            log.info("Annotations before this optimization => ${orderedAnnotationsString}")

                            ArrayList ordAnnArrList = orderedAnnotationsString.split(",")
                            def leastInformativeAnn = ordAnnArrList.remove(ordAnnArrList.size() - 1)

                            def reducedAnnotationsString = ordAnnArrList.toString().replace('[', ' ').replace(']', ' ')
                            log.info("Annotations after eliminating the least informative annotation (${leastInformativeAnn}) => ${reducedAnnotationsString}")
                            return  reducedAnnotationsString
                        }
                }
        }






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





        ann6_ch = GATK_VARIANT_RECALIBRATOR__ANN7.out.annotationsLog
                  .map { logfile ->  eliminateLeastInformativeAnnotation(logfile) }
                  .view()


        GATK_VARIANT_RECALIBRATOR__ANN6(analysisType,
                                    ann6_ch,
                                    select_variants_vcftuple_ch,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )




    emit:

        optimized_vqsr_ch = select_variants_vcftuple_ch
                            .join(GATK_VARIANT_RECALIBRATOR__ANN6.out.recalVcfTuple)
                            .join(GATK_VARIANT_RECALIBRATOR__ANN6.out.tranchesFile)


}
