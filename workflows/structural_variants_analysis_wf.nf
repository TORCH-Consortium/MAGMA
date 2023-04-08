include { TBPROFILER_PROFILE__BAM } from "../modules/tbprofiler/profile__bam.nf" addParams (params.TBPROFILER_PROFILE__BAM)

workflow STRUCTURAL_VARIANTS_ANALYSIS_WF {

    take:
        samtools_bams_ch

    main:

        
        bams_ch = samtools_bams_ch.flatten().collate(3)
         
        bams_ch.dump(tag:"STRUCTURAL_VARIANTS_ANALYSIS_WF")

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        TBPROFILER_PROFILE__BAM(bams_ch, resistanceDb)

}
