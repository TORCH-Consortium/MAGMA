nextflow.enable.dsl = 2

include { GATK_HAPLOTYPE_CALLER } from "../../modules/gatk/haplotype_caller.nf" addParams(params.GATK_HAPLOTYPE_CALLER)

workflow CALL_WF {
    // call_merge
    SAMTOOLS_MERGE

    // call_mark_duplicates
    GATK_MARK_DUPLICATES

    // call_base_recal
    GATK_BASE_RECALIBRATOR

    // call_apply_bqsr
    GATK_APPLY_BQSR
    SAMTOOLS_INDEX

    // call_haplotype_caller
    GATK_HAPLOTYPE_CALLER

    // call_haplotype_caller_minor_variants
    GATK_HAPLOTYPE_CALLER_MINOR_VARIANTS

    // call_ntm
    LOFREQ_CALL

    // call_lofreq
    LOFREQ_INDELQUAL
    SAMTOOLS_INDEX
    LOFREQ_CALL
    LOFREQ_FILTER

    // call_sv
    DELLY_CALL
    BCFTOOLS_VIEW
    GATK_INDEX_FEATURE_FILE
    //TODO: Confirm whether this is correct? (XBS_call#L140)
    //$JAVA -Xmx64G -jar $GATK -V PREFIX.potentialSV.vcf.gz -O PREFIX.potentialSV.DRgenes.vcf.gz -XL $RESOURCE_PATH/DRgenes.list

    // call_stats
    SAMTOOLS_STATS
    GATK_COLLECT_WGS_METRICS
    GATK_FLAG_STAT

}
