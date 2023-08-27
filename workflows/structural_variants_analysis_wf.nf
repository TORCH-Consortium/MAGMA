include { BGZIP } from "../modules/bgzip/bgzip.nf" addParams( params.BGZIP__MINOR_VARIANTS )
include { DELLY_CALL } from "../modules/delly/call.nf" addParams ( params.DELLY_CALL )
include { BWA_MEM as  BWA_MEM__DELLY } from '../modules/bwa/mem.nf' addParams (params.BWA_MEM__DELLY)
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE__DELLY } from "../modules/samtools/merge.nf" addParams ( params.SAMTOOLS_MERGE__DELLY )
include { GATK_MARK_DUPLICATES as  GATK_MARK_DUPLICATES__DELLY } from "../modules/gatk/mark_duplicates.nf" addParams ( params.GATK_MARK_DUPLICATES )
include { GATK_BASE_RECALIBRATOR as GATK_BASE_RECALIBRATOR__DELLY } from "../modules/gatk/base_recalibrator.nf" addParams ( params.GATK_BASE_RECALIBRATOR__DELLY )
include { GATK_APPLY_BQSR as GATK_APPLY_BQSR__DELLY } from "../modules/gatk/apply_bqsr.nf" addParams ( params.GATK_APPLY_BQSR__DELLY )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX__DELLY } from "../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX__DELLY )
include { BCFTOOLS_VIEW__TBP } from "../modules/bcftools/view__tbp.nf" addParams ( params.BCFTOOLS_VIEW__TBP )
include { BCFTOOLS_MERGE as  BCFTOOLS_MERGE__DELLY } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE__DELLY )
include { TBPROFILER_VCF_PROFILE__COHORT as TBPROFILER_VCF_PROFILE__DELLY } from "../modules/tbprofiler/vcf_profile__cohort.nf" addParams (params.TBPROFILER_VCF_PROFILE__DELLY)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__DELLY } from "../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__DELLY)

//include { TBPROFILER_PROFILE__BAM } from "../modules/tbprofiler/profile__bam.nf" addParams (params.TBPROFILER_PROFILE__BAM)


workflow STRUCTURAL_VARIANTS_ANALYSIS_WF {

    take:
        validated_reads_ch
	approved_samples_ch

    main:


        BWA_MEM__DELLY(validated_reads_ch,
                  params.ref_fasta,
                  [params.ref_fasta_dict,
                   params.ref_fasta_amb,
                   params.ref_fasta_ann,
                   params.ref_fasta_bwt,
                   params.ref_fasta_fai,
                   params.ref_fasta_pac,
                   params.ref_fasta_sa])



        normalize_libraries_ch = BWA_MEM__DELLY.out
                                        .map { it -> {
                                                def splittedNameArray = it[0].split("\\.")
                                                def identifier = splittedNameArray[0] + "."  + splittedNameArray[1]

                                                return [identifier, it[1]]
            }
        }
        .groupTuple()
        //.dump(tag: "CALL_WF normalize_libraries_ch : ", pretty: true)

	normalize_filtered_ch = approved_samples_ch.join(normalize_libraries_ch)

        // call_merge
        SAMTOOLS_MERGE__DELLY(normalize_filtered_ch)

        // call_mark_duplicates
        GATK_MARK_DUPLICATES__DELLY(SAMTOOLS_MERGE__DELLY.out)

        if (params.dataset_is_not_contaminated) {
            // call_base_recal
            GATK_BASE_RECALIBRATOR__DELLY(GATK_MARK_DUPLICATES__DELLY.out.bam_tuple,
                                params.dbsnp_vcf,
                                params.ref_fasta,
                                [params.ref_fasta_fai, params.ref_fasta_dict, params.dbsnp_vcf_tbi ] )

            // call_apply_bqsr
            GATK_APPLY_BQSR__DELLY(GATK_BASE_RECALIBRATOR__DELLY.out, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])


            recalibrated_bam_ch = GATK_APPLY_BQSR__DELLY.out

        } else {

            recalibrated_bam_ch = GATK_MARK_DUPLICATES__DELLY.out.bam_tuple
        }


        //recalibrated_bam_ch.dump(tag: "CALL_WF recalibrated_bam_ch: ", pretty:true)

        SAMTOOLS_INDEX__DELLY(recalibrated_bam_ch)


        //----------------------------------------------------------------------------------
        // Infer structural variants
        //
        //NOTE: This inference can not handle a contaminant and MTB allele on the same site.
        //if so the site will be excluded.
        //----------------------------------------------------------------------------------

        DELLY_CALL(SAMTOOLS_INDEX__DELLY.out, params.ref_fasta)

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

	vcfs_string_ch.view() // FIXME I need to view this channel or I get an "Comparison method violates its general contract!" error
	vcfs_file = vcfs_string_ch.collectFile(name: 'structural_variant_vcfs.txt', newLine: true)
        BCFTOOLS_MERGE__DELLY(vcfs_file, vcfs_and_indexes_ch)

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        TBPROFILER_VCF_PROFILE__DELLY(BCFTOOLS_MERGE__DELLY.out, resistanceDb)

	TBPROFILER_COLLATE__DELLY(params.vcf_name, TBPROFILER_VCF_PROFILE__DELLY.out, resistanceDb)
	
    emit:
	structural_variants_results_ch = TBPROFILER_COLLATE__DELLY.out.per_sample_results
}
