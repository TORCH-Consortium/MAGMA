include { TBPROFILER_PROFILER__BAM } from "../modules/tbprofiler/profile__bam.nf" addParams (params.TBPROFILER_PROFILER__BAM)

workflow STRUCTURAL_VARIANTS_ANALYSIS_WF {

    take:
        samtools_bams_ch

    main:

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []
        samtools_bams_ch.dump(tag:"STRUCTURAL_VARIANTS_ANALYSIS_WF")

        TBPROFILER_PROFILER__BAM(samtools_bams_ch, resistanceDb)

}
