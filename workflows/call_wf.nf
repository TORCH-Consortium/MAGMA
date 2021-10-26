include { SAMTOOLS_MERGE } from "../modules/samtools/merge.nf" addParams ( params.SAMTOOLS_MERGE )
// include { GATK_MARK_DUPLICATES } from "../modules/gatk/mark_duplicates.nf" addParams ( params.GATK_MARK_DUPLICATES )
// include { GATK_BASE_RECALIBRATOR } from "../modules/gatk/base_recalibrator.nf" addParams ( params.GATK_BASE_RECALIBRATOR )
// include { GATK_APPLY_BQSR } from "../modules/gatk/apply_bqsr.nf" addParams ( params.GATK_APPLY_BQSR )
// include { SAMTOOLS_INDEX } from "../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX )
// include { GATK_HAPLOTYPE_CALLER } from "../modules/gatk/haplotype_caller.nf" addParams ( params.GATK_HAPLOTYPE_CALLER )
// include { GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS } from "../modules/gatk/haplotype_caller__minor_variants.nf" addParams ( params.GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS )
// include { LOFREQ_CALL__NTM } from "../modules/lofreq/call__ntm.nf" addParams ( params.LOFREQ_CALL__NTM )
// include { LOFREQ_INDELQUAL } from "../modules/lofreq/indelqual.nf" addParams ( params.LOFREQ_INDELQUAL )
// include { SAMTOOLS_INDEX as SAMTOOLS_INDEX__LOFREQ } from "../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX__LOFREQ )
// include { LOFREQ_CALL } from "../modules/lofreq/call.nf" addParams ( params.LOFREQ_CALL )
// include { LOFREQ_FILTER } from "../modules/lofreq/filter.nf" addParams ( params.LOFREQ_FILTER )
// include { DELLY_CALL } from "../modules/delly/call.nf" addParams ( params.DELLY_CALL )
// include { BCFTOOLS_VIEW } from "../modules/bcftools/view.nf" addParams ( params.BCFTOOLS_VIEW )
// include { GATK_INDEX_FEATURE_FILE } from "../modules/gatk/index_feature_file.nf" addParams ( params.GATK_INDEX_FEATURE_FILE )
// include { SAMTOOLS_STATS } from "../modules/samtools/stats.nf" addParams ( params.SAMTOOLS_STATS )
// include { GATK_COLLECT_WGS_METRICS } from "../modules/gatk/collect_wgs_metrics.nf" addParams ( params.GATK_COLLECT_WGS_METRICS )
// include { GATK_FLAG_STAT } from "../modules/gatk/flag_stat.nf" addParams ( params.GATK_FLAG_STAT )
// include { UTILS_SAMPLE_STATS } from "../modules/utils/sample_stats.nf" addParams ( params.UTILS_SAMPLE_STATS )
// include { UTILS_COHORT_STATS } from "../modules/utils/cohort_stats.nf" addParams ( params.UTILS_COHORT_STATS )





workflow CALL_WF {
    take:
        sorted_reads_ch


    main:

        // call_merge
        //NOTE: The output of this seems to overwrite the output of XBS_map#L31
        //TODO: Confirm whether this is really necessary as the results are the same and we don't have multiple BAM files
        SAMTOOLS_MERGE(sorted_reads_ch)

        // call_mark_duplicates
        GATK_MARK_DUPLICATES(SAMTOOLS_MERGE.out)

    /*
        if (params.dataset_is_not_contaminated) {
            // call_base_recal
            GATK_BASE_RECALIBRATOR(GATK_MARK_DUPLICATES.out, params.dbsnp_vcf, params.ref_fasta)

            // call_apply_bqsr
            GATK_APPLY_BQSR(GATK_BASE_RECALIBRATOR.out, params.ref_fasta)


            recalibrated_bam_ch = GATK_APPLY_BQSR.out

        } else {

            recalibrated_bam_ch = GATK_MARK_DUPLICATES.out
        }

        SAMTOOLS_INDEX(recalibrated_bam_ch)

        //----------------------------------------------------------------------------------
        // Call Variants for follow up calling
        //----------------------------------------------------------------------------------



        // call_haplotype_caller
        GATK_HAPLOTYPE_CALLER(SAMTOOLS_INDEX.out, params.ref_fasta)

        if (params.compute_minor_variants) {
            // call_haplotype_caller_minor_variants
            GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS(SAMTOOLS_INDEX.out, params.ref_fasta)
        }

        //----------------------------------------------------------------------------------
        // Infer potential NTM contamination
        //----------------------------------------------------------------------------------

        // call_ntm
        LOFREQ_CALL__NTM(GATK_HAPLOTYPE_CALLER.out, params.ref_fasta)


        //----------------------------------------------------------------------------------
        // Infer minor variants with LoFreq
        //----------------------------------------------------------------------------------

        // call_lofreq
        LOFREQ_INDELQUAL(GATK_HAPLOTYPE_CALLER.out, params.ref_fasta)
        SAMTOOLS_INDEX__LOFREQ(LOFREQ_INDELQUAL.out,params.ref_fasta)
        LOFREQ_CALL(SAMTOOLS_INDEX.out)
        LOFREQ_FILTER(LOFREQ_CALL.out)

        //----------------------------------------------------------------------------------
        // Infer structural variants
        //
        //NOTE: This inference can not handle a contaminant and MTB allele on the same site, if so the site will be excluded.
        //----------------------------------------------------------------------------------

        // call_sv
        //TODO: Should SV-calling be optional?

        DELLY_CALL(recalibrated_bam_ch)
        BCFTOOLS_VIEW(DELLY_CALL.out)
        GATK_INDEX_FEATURE_FILE(BCFTOOLS_VIEW.out)

        //Enable this once a proper file with DR genes has been made:
        GATK_SELECT_VARIANTS__INTERVALS(GATK_INDEX_FEATURE_FILE.out, params.drgenes_list)


        //----------------------------------------------------------------------------------
        // STATS
        //----------------------------------------------------------------------------------


        // call_stats
        SAMTOOLS_STATS(recalibrated_bam_ch, params.ref_fasta)
        GATK_COLLECT_WGS_METRICS(recalibrated_bam_ch, params.ref_fasta)
        GATK_FLAG_STAT(recalibrated_bam_ch, params.ref_fasta)


        sample_stats_ch = (SAMTOOLS_STATS.out)
            .join(GATK_COLLECT_WGS_METRICS.out)
            .join(GATK_FLAG_STAT.out)
            .join(QUANTTB_QUANT.out)
            .join(LOFREQ_CALL__NTM.out)


        UTILS_SAMPLE_STATS(sample_stats_ch)

        UTILS_COHORT_STATS(UTILS_SAMPLE_STATS.out.collect())

    emit:
        cohort_stats_tsv = UTILS_COHORT_STATS.out
*/
}
