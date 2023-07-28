include { BGZIP } from "../modules/bgzip/bgzip.nf" addParams( params.BGZIP__MINOR_VARIANTS )
include { DELLY_CALL } from "../modules/delly/call.nf" addParams ( params.DELLY_CALL )
include { BCFTOOLS_VIEW__GATK} from "../modules/bcftools/view__gatk.nf" addParams ( params.BCFTOOLS_VIEW__GATK )
include { BCFTOOLS_VIEW__TBP } from "../modules/bcftools/view__tbp.nf" addParams ( params.BCFTOOLS_VIEW__TBP )
include { BCFTOOLS_MERGE } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE )
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE__SV } from "../modules/gatk/index_feature_file.nf" addParams ( params.GATK_INDEX_FEATURE_FILE__SV )
include { GATK_SELECT_VARIANTS__INCLUSION } from "../modules/gatk/select_variants__intervals.nf" addParams ( params.GATK_SELECT_VARIANTS__INCLUSION )
include { TBPROFILER_PROFILE__BAM } from "../modules/tbprofiler/profile__bam.nf" addParams (params.TBPROFILER_PROFILE__BAM)

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

        DELLY_CALL(samtools_bams_ch, params.ref_fasta)

        BCFTOOLS_VIEW__GATK(DELLY_CALL.out)
        GATK_INDEX_FEATURE_FILE__SV(BCFTOOLS_VIEW__GATK.out, 'potentialSV')
        GATK_SELECT_VARIANTS__INCLUSION(GATK_INDEX_FEATURE_FILE__SV.out.sample_vcf_tuple, params.drgenes_list)

        BCFTOOLS_VIEW__TBP(DELLY_CALL.out)

//FIXME save the string to an intermediate file
        vcfs_string_ch = BCFTOOLS_VIEW__TBP.out
                                .collect()
                                .dump(tag:'MINOR_VARIANT_WF: vcfs_string_ch', pretty: true)

                                /*
                                .flatten()
                                .filter { it.extension  == "gz" }
                                .map { it -> it.name }
                                .reduce { a, b -> "$a $b " }
                                .dump(tag:'MINOR_VARIANT_WF: vcfs_string_ch', pretty: true)
                                */


        //BCFTOOLS_MERGE(vcfs_string_ch, reformatted_lofreq_vcfs_tuple_ch)

        // merge_call_resistance_lofreq
        //BGZIP(BCFTOOLS_MERGE.out) 

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        bams_ch = samtools_bams_ch
                        .collect()
                        .flatten()
                        .collate(3)
                        .dump(tag:"STRUCTURAL_VARIANTS_ANALYSIS_WF")


        TBPROFILER_PROFILE__BAM(bams_ch, resistanceDb)

}
