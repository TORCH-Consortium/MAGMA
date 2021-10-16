
workflow RESISTANCE_ANALYSIS {

    // merge_call_resistance
    TBPROFILER_VCF_PROFILE__COHORT
    TBPROFILER_COLLATE

    // merge_call_resistance_lofreq
    BGZIP()
    TBPROFILER_VCF_PROFILE__SAMPLE

    TBPROFILER_COLLATE()

}
