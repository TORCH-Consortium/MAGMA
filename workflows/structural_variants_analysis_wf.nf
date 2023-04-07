include { TBPROFILER_PROFILE as TBPROFILER_PROFILER__BAM } from "../modules/tbprofiler/profile.nf" addParams (params.TBPROFILER_PROFILE)

workflow STRUCTURAL_VARIANTS_ANALYSIS_WF {

    take:
        samtools_bams_ch

    main:

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        TBPROFILER_PROFILER__BAM(samtools_bams_ch, resistanceDb)

}
