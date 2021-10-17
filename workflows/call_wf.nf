
include { SAMTOOLS_MERGE } from "../../modules/samtools/merge.nf" addParams ( params.SAMTOOLS_MERGE )
include { GATK_MARK_DUPLICATES } from "../../modules/gatk/mark_duplicates.nf" addParams ( params.GATK_MARK_DUPLICATES )
include { GATK_BASE_RECALIBRATOR } from "../../modules/gatk/base_recalibrator.nf" addParams ( params.GATK_BASE_RECALIBRATOR )
include { GATK_APPLY_BQSR } from "../../modules/gatk/apply_bqsr.nf" addParams ( params.GATK_APPLY_BQSR )
include { SAMTOOLS_INDEX } from "../../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX )
include { GATK_HAPLOTYPE_CALLER } from "../../modules/gatk/haplotype_caller.nf" addParams ( params.GATK_HAPLOTYPE_CALLER )
include { GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS } from "../../modules/gatk/haplotype_caller__minor_variants.nf" addParams ( params.GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS )
include { LOFREQ_CALL__NTM } from "../../modules/lofreq/call__ntm.nf" addParams ( params.LOFREQ_CALL__NTM )
include { LOFREQ_INDELQUAL } from "../../modules/lofreq/indelqual.nf" addParams ( params.LOFREQ_INDELQUAL )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX__LOFREQ } from "../../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX__LOFREQ )
include { LOFREQ_CALL } from "../../modules/lofreq/call.nf" addParams ( params.LOFREQ_CALL )
include { LOFREQ_FILTER } from "../../modules/lofreq/fileter.nf" addParams ( params.LOFREQ_FILTER )
include { DELLY_CALL } from "../../modules/delly/call.nf" addParams ( params.DELLY_CALL )
include { BCFTOOLS_VIEW } from "../../modules/bcftools/view.nf" addParams ( params.BCFTOOLS_VIEW )
include { GATK_INDEX_FEATURE_FILE } from "../../modules/gatk/index_feature_file.nf" addParams ( params.GATK_INDEX_FEATURE_FILE )
include { SAMTOOLS_STATS } from "../../modules/samtools/stats.nf" addParams ( params.SAMTOOLS_STATS )
include { GATK_COLLECT_WGS_METRICS } from "../../modules/gatk/collect_wgs_metrics.nf" addParams ( params.GATK_COLLECT_WGS_METRICS )
include { GATK_FLAG_STAT } from "../../modules/gatk/flag_stat.nf" addParams ( params.GATK_FLAG_STAT )
include { UTILS_SAMPLE_STATS } from "../../modules/utils/sample_stats.nf" addParams ( params.UTILS_SAMPLE_STATS )
include { UTILS_COHORT_STATS } from "../../modules/utils/cohort_stats.nf" addParams ( params.UTILS_COHORT_STATS )





workflow CALL_WF {
    take:


    main:

        // call_merge
        //NOTE: The output of this seems to overwrite the output of XBS_map#L31
        SAMTOOLS_MERGE()

        // call_mark_duplicates
        GATK_MARK_DUPLICATES(SAMTOOLS_MERGE.out)

        // call_base_recal
        //NOTE: (XBS_call#L48) Disabled by default - Enable for non-contaminated datasets only.
        //GATK_BASE_RECALIBRATOR(GATK_MARK_DUPLICATES.out, params.dbsnp, params.ref_fasta)

        // call_apply_bqsr
        //NOTE: (XBS_call#L63) Seems to depend up on the BASE_RECALIBRATOR - commented out by default.
        // GATK_APPLY_BQSR(GATK_BASE_RECALIBRATOR.out, params.ref_fasta)

        //NOTE: If APPLY_BQSR is NOT used, then this directly connects to MARK_DUPLICATES
        SAMTOOLS_INDEX(GATK_MARK_DUPLICATES.out)

        // call_haplotype_caller
        GATK_HAPLOTYPE_CALLER(SAMTOOLS_INDEX.out, params.ref_fasta)

        // call_haplotype_caller_minor_variants
        GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS(SAMTOOLS_INDEX.out, params.ref_fasta)

        // call_ntm
        LOFREQ_CALL__NTM(GATK_HAPLOTYPE_CALLER.out, params.ref_fasta)

        // call_lofreq
        LOFREQ_INDELQUAL(GATK_HAPLOTYPE_CALLER.out, params.ref_fasta)
        SAMTOOLS_INDEX__LOFREQ(LOFREQ_INDELQUAL.out,params.ref_fasta)
        LOFREQ_CALL(SAMTOOLS_INDEX.out)
        LOFREQ_FILTER(LOFREQ_CALL.out)

        // call_sv
        DELLY_CALL(GATK_MARK_DUPLICATES.out)
        BCFTOOLS_VIEW(DELLY_CALL.out)
        GATK_INDEX_FEATURE_FILE(BCFTOOLS_VIEW)

        //Enable this once a proper file with DR genes has been made:
        //FIXME add inputs
        GATK_SELECT_VARIANTS__INTERVALS()

        // FIXME add the inputs
        // call_stats
        SAMTOOLS_STATS
        GATK_COLLECT_WGS_METRICS
        GATK_FLAG_STAT


        //TODO: Make sure that all of the stat files belong to the same sample (join operator)
        UTILS_SAMPLE_STATS

        UTILS_COHORT_STATS(UTILS_SAMPLE_STATS.out.collect())

}
