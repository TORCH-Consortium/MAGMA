include { TBPROFILER_VCF_PROFILE__COHORT } from "../../modules/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__COHORT)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__COHORT } from "../../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__COHORT)

workflow RESISTANCE_ANALYSIS {
    take:
        merged_vcf_ch
        reformatted_lofreq_vcf_ch

    main:

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        // merge_call_resistance
        TBPROFILER_VCF_PROFILE__COHORT(merged_vcf_ch, resistanceDb)
        TBPROFILER_COLLATE__COHORT(params.vcf_name, TBPROFILER_VCF_PROFILE__COHORT.out, resistanceDb)

        // merge_call_resistance_lofreq
        // NOTE: Moved to the MINOR_COHORT_VARIANT_WF
        /* BGZIP(lofreq_vcf_ch) */
        /* TBPROFILER_VCF_PROFILE__LOFREQ(BGZIP.out, resistanceDb) */
        /* TBPROFILER_COLLATE__LOFREQ(params.vcf_name, */
        /*                           TBPROFILER_VCF_PROFILE__LOFREQ.out.resistance_json.collect(), */
        /*                           resistanceDb) */
}
