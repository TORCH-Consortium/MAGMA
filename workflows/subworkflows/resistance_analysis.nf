include { TBPROFILER_VCF_PROFILE__COHORT } from "../../modules/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__COHORT)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__COHORT } from "../../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__COHORT)


workflow RESISTANCE_ANALYSIS {
    take:
        merged_vcf_ch

    main:

        def database = params.resistance_db ? params.resistance_db : []

        // merge_call_resistance
        TBPROFILER_VCF_PROFILE__COHORT(merged_vcf_ch, database)

        TBPROFILER_COLLATE__COHORT(TBPROFILER_VCF_PROFILE__COHORT.out, database)

        /*
        // merge_call_resistance_lofreq
        BGZIP()
        TBPROFILER_VCF_PROFILE__LOFREQ

        TBPROFILER_COLLATE__LOFREQ()

        */
}
