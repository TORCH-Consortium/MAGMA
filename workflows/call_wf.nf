include { BGZIP as BGZIP__LOFREQ } from "../modules/bgzip/bgzip.nf" addParams( params.BGZIP__LOFREQ )
include { SAMTOOLS_MERGE } from "../modules/samtools/merge.nf" addParams ( params.SAMTOOLS_MERGE )
include { GATK_MARK_DUPLICATES } from "../modules/gatk/mark_duplicates.nf" addParams ( params.GATK_MARK_DUPLICATES )
include { GATK_BASE_RECALIBRATOR } from "../modules/gatk/base_recalibrator.nf" addParams ( params.GATK_BASE_RECALIBRATOR )
include { GATK_APPLY_BQSR } from "../modules/gatk/apply_bqsr.nf" addParams ( params.GATK_APPLY_BQSR )
include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller.nf" addParams ( params.GATK_HAPLOTYPE_CALLER )
include { GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS } from "../modules/gatk/haplotype_caller__minor_variants.nf" addParams ( params.GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS )
include { LOFREQ_CALL__NTM } from "../modules/lofreq/call__ntm.nf" addParams ( params.LOFREQ_CALL__NTM )
include { LOFREQ_INDELQUAL } from "../modules/lofreq/indelqual.nf" addParams ( params.LOFREQ_INDELQUAL )
include { SAMTOOLS_INDEX } from "../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX )
include { SAMTOOLS_INDEX__LOFREQ } from "../modules/samtools/index__lofreq.nf" addParams ( params.SAMTOOLS_INDEX__LOFREQ )
include { LOFREQ_CALL } from "../modules/lofreq/call.nf" addParams ( params.LOFREQ_CALL )
include { LOFREQ_FILTER } from "../modules/lofreq/filter.nf" addParams ( params.LOFREQ_FILTER )
include { SAMTOOLS_STATS } from "../modules/samtools/stats.nf" addParams ( params.SAMTOOLS_STATS )
include { GATK_COLLECT_WGS_METRICS } from "../modules/gatk/collect_wgs_metrics.nf" addParams ( params.GATK_COLLECT_WGS_METRICS )
include { GATK_FLAG_STAT } from "../modules/gatk/flag_stat.nf" addParams ( params.GATK_FLAG_STAT )
include { UTILS_SAMPLE_STATS } from "../modules/utils/sample_stats.nf" addParams ( params.UTILS_SAMPLE_STATS )
include { UTILS_COHORT_STATS } from "../modules/utils/cohort_stats.nf" addParams ( params.UTILS_COHORT_STATS )
include { UTILS_REFORMAT_LOFREQ } from "../modules/utils/reformat_lofreq.nf" addParams ( params.UTILS_REFORMAT_LOFREQ )
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE__LOFREQ } from "../modules/gatk/index_feature_file.nf" addParams ( params.GATK_INDEX_FEATURE_FILE__LOFREQ )



workflow CALL_WF {
    take:
        bam_sorted_reads_ch


    main:

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

        if (!params.skip_base_recalibration) {
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
        // Call Variants for follow up calling
        //----------------------------------------------------------------------------------



        // call_haplotype_caller
        GATK_HAPLOTYPE_CALLER(SAMTOOLS_INDEX.out,
                          params.ref_fasta,
                          [params.ref_fasta_fai, params.ref_fasta_dict])

        // call_haplotype_caller_minor_variants
        if (!params.skip_minor_variants_gatk) {
            GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS(SAMTOOLS_INDEX.out,
                                                params.ref_fasta,
                                                [params.ref_fasta_fai, params.ref_fasta_dict])
        }

        //----------------------------------------------------------------------------------
        // Infer potential NTM contamination
        //----------------------------------------------------------------------------------


        // call_ntm
        LOFREQ_CALL__NTM(SAMTOOLS_INDEX.out,
                         params.ref_fasta,
                         [params.ref_fasta_fai])

        //----------------------------------------------------------------------------------
        // Infer minor variants with LoFreq
        //----------------------------------------------------------------------------------

        // call_lofreq
        LOFREQ_INDELQUAL(recalibrated_bam_ch, params.ref_fasta)
        SAMTOOLS_INDEX__LOFREQ(LOFREQ_INDELQUAL.out)
        LOFREQ_CALL(SAMTOOLS_INDEX__LOFREQ.out, params.ref_fasta, [params.ref_fasta_fai])
        LOFREQ_FILTER(LOFREQ_CALL.out, params.ref_fasta)

        //----------------------------------------------------------------------------------
        // Reformat individual results from LOFREQ
        //----------------------------------------------------------------------------------

        UTILS_REFORMAT_LOFREQ(LOFREQ_CALL.out)

        BGZIP__LOFREQ(UTILS_REFORMAT_LOFREQ.out)

        GATK_INDEX_FEATURE_FILE__LOFREQ(BGZIP__LOFREQ.out, 'lofreq.reformat.corrected')

        //----------------------------------------------------------------------------------
        // STATS
        //----------------------------------------------------------------------------------


        // call_stats
        SAMTOOLS_STATS(recalibrated_bam_ch, params.ref_fasta)
        GATK_COLLECT_WGS_METRICS(recalibrated_bam_ch, params.ref_fasta)
        GATK_FLAG_STAT(recalibrated_bam_ch, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])


        sample_stats_ch = (SAMTOOLS_STATS.out)
            .join(GATK_COLLECT_WGS_METRICS.out)
            .join(GATK_FLAG_STAT.out)
            .join(LOFREQ_CALL__NTM.out)
            //.dump(tag: "CALL_WF sample_stats_ch : ", pretty: true)


        UTILS_SAMPLE_STATS(sample_stats_ch)

        UTILS_COHORT_STATS(UTILS_SAMPLE_STATS.out.collect())

    emit:
        cohort_stats_tsv = UTILS_COHORT_STATS.out
        gvcf_ch = GATK_HAPLOTYPE_CALLER.out.gvcf_ch.collect()
        reformatted_lofreq_vcfs_tuple_ch = GATK_INDEX_FEATURE_FILE__LOFREQ.out.vcf_tuple.collect(sort:true)
        bgzip_ch = BGZIP__LOFREQ.out.collect() 
        samtools_bam_ch = SAMTOOLS_INDEX.out 
}
