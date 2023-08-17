include { BGZIP } from "../modules/bgzip/bgzip.nf" addParams( params.BGZIP__MINOR_VARIANTS )
include { DELLY_CALL } from "../modules/delly/call.nf" addParams ( params.DELLY_CALL )
include { BCFTOOLS_VIEW__GATK} from "../modules/bcftools/view__gatk.nf" addParams ( params.BCFTOOLS_VIEW__GATK )
include { BCFTOOLS_VIEW__TBP } from "../modules/bcftools/view__tbp.nf" addParams ( params.BCFTOOLS_VIEW__TBP )
include { BCFTOOLS_MERGE as  BCFTOOLS_MERGE__DELLY } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE__DELLY )
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE__SV } from "../modules/gatk/index_feature_file.nf" addParams ( params.GATK_INDEX_FEATURE_FILE__SV )
include { GATK_SELECT_VARIANTS__INCLUSION } from "../modules/gatk/select_variants__intervals.nf" addParams ( params.GATK_SELECT_VARIANTS__INCLUSION )
include { TBPROFILER_VCF_PROFILE__COHORT as TBPROFILER_VCF_PROFILE__DELLY } from "../modules/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__DELLY)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__DELLY } from "../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__DELLY)

//include { TBPROFILER_PROFILE__BAM } from "../modules/tbprofiler/profile__bam.nf" addParams (params.TBPROFILER_PROFILE__BAM)

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

        BCFTOOLS_VIEW__TBP(DELLY_CALL.out)

	//FIXME save the string to an intermediate file

        vcfs_and_indexes_ch = BCFTOOLS_VIEW__TBP.out
                                .collect(sort: true)
                                .flatten()
                                .filter { it.class.name  != "java.lang.String" }
                                .collect(sort: true)
                                //.view{ it }

        vcfs_string_ch = BCFTOOLS_VIEW__TBP.out
                                .collect(sort: true)
                                .flatten()
                                .filter { it.class.name  != "java.lang.String" }
                                .filter { it.extension  == "gz" }
                                .map { it -> it.name }
                                .reduce { a, b -> "$a $b " }
                                //.view { it }
                                //.dump(tag:'MINOR_VARIANT_WF: vcfs_string_ch', pretty: true)


        BCFTOOLS_MERGE__DELLY(vcfs_string_ch, vcfs_and_indexes_ch)

	def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        TBPROFILER_VCF_PROFILE__DELLY(BCFTOOLS_MERGE__DELLY.out, resistanceDb)

	TBPROFILER_COLLATE__DELLY(params.vcf_name, TBPROFILER_VCF_PROFILE__DELLY.out, resistanceDb)
	
    emit:
	structural_variants_results_ch = TBPROFILER_COLLATE__DELLY.out.per_sample_results
}
