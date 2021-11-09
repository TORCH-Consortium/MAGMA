include { TBPROFILER_VCF_PROFILE__COHORT } from "../../modules/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__COHORT)


workflow RESISTANCE_ANALYSIS {
    take:
        merged_vcf_ch

    main:

        database = params.resistance_db ? params.resistance_db : []

        // merge_call_resistance
        TBPROFILER_VCF_PROFILE__COHORT(merged_vcf_ch, database)

        // TBPROFILER_COLLATE

        /*
        // merge_call_resistance_lofreq
        BGZIP()
        TBPROFILER_VCF_PROFILE__SAMPLE

        TBPROFILER_COLLATE()

        */
}
