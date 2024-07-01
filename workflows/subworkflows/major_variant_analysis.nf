include { TBPROFILER_VCF_PROFILE__COHORT } from "../../modules/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__COHORT)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__COHORT } from "../../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__COHORT)

workflow MAJOR_VARIANT_ANALYSIS {
    take:
        merged_vcf_ch
        reformatted_lofreq_vcf_ch

    main:

        def resistanceDb =  []

        // merge_call_resistance
        TBPROFILER_VCF_PROFILE__COHORT(merged_vcf_ch, resistanceDb)
        TBPROFILER_COLLATE__COHORT(params.vcf_name, TBPROFILER_VCF_PROFILE__COHORT.out, resistanceDb)

        emit:
            major_variants_results_ch =  TBPROFILER_COLLATE__COHORT.out.per_sample_results
}
