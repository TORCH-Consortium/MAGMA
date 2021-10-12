nextflow.enable.dsl = 2

include { GATK_HAPLOTYPE_CALLER } from "../../modules/gatk/haplotype_caller.nf" addParams(params.GATK_HAPLOTYPE_CALLER)

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

        // NOTE: If APPLY_BQSR is NOT used, then this directly connects to MARK_DUPLICATES
        SAMTOOLS_INDEX(GATK_MARK_DUPLICATES.out)

        // call_haplotype_caller
        GATK_HAPLOTYPE_CALLER(SAMTOOLS_INDEX.out, params.ref_fasta)

        //TODO
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
