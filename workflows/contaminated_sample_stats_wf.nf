include { GATK_MARK_DUPLICATES } from "../modules/gatk/mark_duplicates.nf" addParams ( params.GATK_MARK_DUPLICATES__REJECTED )
include { SAMTOOLS_INDEX } from "../modules/samtools/index.nf" addParams ( params.SAMTOOLS_INDEX__REJECTED )
include { SAMTOOLS_STATS } from "../modules/samtools/stats.nf" addParams ( params.SAMTOOLS_STATS__REJECTED )
include { GATK_COLLECT_WGS_METRICS } from "../modules/gatk/collect_wgs_metrics.nf" addParams ( params.GATK_COLLECT_WGS_METRICS__REJECTED )
include { GATK_FLAG_STAT } from "../modules/gatk/flag_stat.nf" addParams ( params.GATK_FLAG_STAT__REJECTED )
include { UTILS_SAMPLE_STATS } from "../modules/utils/sample_stats.nf" addParams ( params.UTILS_SAMPLE_STATS__REJECTED )
include { UTILS_COHORT_STATS } from "../modules/utils/cohort_stats.nf" addParams ( params.UTILS_COHORT_STATS__REJECTED )


workflow CONTAMINATED_SAMPLE_STATS_WF {
    take:
        rejected_sorted_reads_ch

    main:

        GATK_MARK_DUPLICATES( rejected_sorted_reads_ch )

        recalibrated_bam_ch = GATK_MARK_DUPLICATES.out.bam_tuple

        SAMTOOLS_INDEX(recalibrated_bam_ch)

        SAMTOOLS_STATS(recalibrated_bam_ch, params.ref_fasta)
        GATK_COLLECT_WGS_METRICS(recalibrated_bam_ch, params.ref_fasta)
        GATK_FLAG_STAT(recalibrated_bam_ch, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])


    /*
        sample_stats_ch = ( SAMTOOLS_STATS.out )
            .join( GATK_COLLECT_WGS_METRICS.out )
            .join( GATK_FLAG_STAT.out )


        UTILS_SAMPLE_STATS( sample_stats_ch )
        UTILS_COHORT_STATS( UTILS_SAMPLE_STATS.out.collect() )

    */
}
