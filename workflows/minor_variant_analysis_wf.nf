workflow MINOR_VARIANT_ANALYSIS_WF {

/* FIXME */
    take:
        reformatted_lofreq_vcf_ch

    main:

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

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

/* FIXME */
    emit:
        UTILS_MULTIPLE_INFECTION_FILTER

}
