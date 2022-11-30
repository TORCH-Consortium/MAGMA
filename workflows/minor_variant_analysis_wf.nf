include { BCFTOOLS_MERGE } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE )


workflow MINOR_VARIANT_ANALYSIS_WF {

    take:
        reformatted_lofreq_vcfs_tuple_ch

    main:

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        vcfs_string_ch = reformatted_lofreq_vcfs_tuple_ch
                                .flatten()
                                .filter { it.extension  == "gz" }
                                .map { it -> it.name }
                                .reduce { a, b -> "$a $b " }
                                /* .view {"\n\n XBS-NF-LOG vcfs_filenames_ch : $it \n\n"} */

        BCFTOOLS_MERGE(vcfs_string_ch, reformatted_lofreq_vcfs_tuple_ch)

/*
        // merge_call_resistance_lofreq
        BGZIP(lofreq_vcf_ch) 

        BGZIP_COHORT_FILE

        //TBPROFILER minor variants
        TBPROFILER_VCF_PROFILE__LOFREQ(BGZIP_COHORT_FILE.out, resistanceDb)

        //TBPROFILER major variants
        TBPROFILER_COLLATE__LOFREQ(params.vcf_name,
                                  TBPROFILER_VCF_PROFILE__LOFREQ.out.resistance_json.collect(),
                                  resistanceDb)

        //TBPROFILER major variants
        UTILS_MULTIPLE_INFECTION_FILTER

    emit:
        UTILS_MULTIPLE_INFECTION_FILTER

*/
}
