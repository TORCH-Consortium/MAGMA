include { TBPROFILER_PROFILE__BAM } from "../modules/tbprofiler/profile__bam.nf" addParams (params.TBPROFILER_PROFILE__BAM)
include { DELLY_CALL } from "../modules/delly/call.nf" addParams ( params.DELLY_CALL )
include { BCFTOOLS_VIEW } from "../modules/bcftools/view.nf" addParams ( params.BCFTOOLS_VIEW )
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE__SV } from "../modules/gatk/index_feature_file.nf" addParams ( params.GATK_INDEX_FEATURE_FILE__SV )
include { GATK_SELECT_VARIANTS__INCLUSION } from "../modules/gatk/select_variants__intervals.nf" addParams ( params.GATK_SELECT_VARIANTS__INCLUSION )

workflow STRUCTURAL_VARIANTS_ANALYSIS_WF {

    take:
        samtools_bams_ch

    main:

        //----------------------------------------------------------------------------------
        // Infer structural variants
        //
        //NOTE: This inference can not handle a contaminant and MTB allele on the same site.
        //if so the site will be excluded.
        //----------------------------------------------------------------------------------


        // call_sv
        DELLY_CALL(samtools_bams_ch, params.ref_fasta)
        BCFTOOLS_VIEW(DELLY_CALL.out)
        GATK_INDEX_FEATURE_FILE__SV(BCFTOOLS_VIEW.out, 'potentialSV')
        GATK_SELECT_VARIANTS__INCLUSION(GATK_INDEX_FEATURE_FILE__SV.out.sample_vcf_tuple, params.drgenes_list)

        bams_ch = samtools_bams_ch
                        .collect()
                        .flatten()
                        .collate(3)
                        .dump(tag:"STRUCTURAL_VARIANTS_ANALYSIS_WF")

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        TBPROFILER_PROFILE__BAM(bams_ch, resistanceDb)

}
