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


//NOTE: -k 19 for BWA_MEM and downstream SAMTOOLS

workflow STRUCTURAL_VARIANTS_ANALYSIS_WF {

    take:
        validated_reads_ch

    main:


        BWA_MEM__DELLY(samples_ch,
                  params.ref_fasta,
                  [params.ref_fasta_dict,
                   params.ref_fasta_amb,
                   params.ref_fasta_ann,
                   params.ref_fasta_bwt,
                   params.ref_fasta_fai,
                   params.ref_fasta_pac,
                   params.ref_fasta_sa])



        normalize_libraries_ch = bam_sorted_reads_ch
                                        .map { it -> {
                                                def splittedNameArray = it[0].split("\\.")
                                                def identifier = splittedNameArray[0] + "."  + splittedNameArray[1]

                                                return [identifier, it[1]]
            }
        }
        .groupTuple()
        //.dump(tag: "CALL_WF normalize_libraries_ch : ", pretty: true)


        // call_merge
        SAMTOOLS_MERGE(normalize_libraries_ch)

        // call_mark_duplicates
        GATK_MARK_DUPLICATES(SAMTOOLS_MERGE.out)

        if (params.dataset_is_not_contaminated) {
            // call_base_recal
            GATK_BASE_RECALIBRATOR(GATK_MARK_DUPLICATES.out.bam_tuple,
                                params.dbsnp_vcf,
                                params.ref_fasta,
                                [params.ref_fasta_fai, params.ref_fasta_dict, params.dbsnp_vcf_tbi ] )

            // call_apply_bqsr
            GATK_APPLY_BQSR(GATK_BASE_RECALIBRATOR.out, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])


            recalibrated_bam_ch = GATK_APPLY_BQSR.out

        } else {

            recalibrated_bam_ch = GATK_MARK_DUPLICATES.out.bam_tuple
        }


        //recalibrated_bam_ch.dump(tag: "CALL_WF recalibrated_bam_ch: ", pretty:true)

        SAMTOOLS_INDEX(recalibrated_bam_ch)


        //----------------------------------------------------------------------------------
        // Infer structural variants
        //
        //NOTE: This inference can not handle a contaminant and MTB allele on the same site.
        //if so the site will be excluded.
        //----------------------------------------------------------------------------------

        DELLY_CALL(SAMTOOLS_INDEX.out, params.ref_fasta)

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


}
