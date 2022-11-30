include { BCFTOOLS_MERGE } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE )


workflow MINOR_VARIANT_ANALYSIS_WF {

    take:
        reformatted_lofreq_vcfs_tuple_ch

    main:

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        //Filter out only the VCF files from the channel
        vcfs_filenames_ch = reformatted_lofreq_vcfs_tuple_ch
                                .view {"\n\n XBS-NF-LOG vcfs_filenames_ch : ${it.split("\.")[-1]} \n\n"}
                                /* .filter { it.class  == sun.nio.fs.UnixPath } */

/*
        BCFTOOLS_MERGE

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
