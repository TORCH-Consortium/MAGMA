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

        //NOTE: If APPLY_BQSR is NOT used, then this directly connects to MARK_DUPLICATES
        SAMTOOLS_INDEX(GATK_MARK_DUPLICATES.out)

        // call_haplotype_caller
        GATK_HAPLOTYPE_CALLER(SAMTOOLS_INDEX.out, params.ref_fasta)

        // call_haplotype_caller_minor_variants
        GATK_HAPLOTYPE_CALLER__MINOR_VARIANTS(SAMTOOLS_INDEX.out, params.ref_fasta)

        // call_ntm
        LOFREQ_CALL_NTM(GATK_HAPLOTYPE_CALLER.out, params.ref_fasta)

        // call_lofreq
        LOFREQ_INDELQUAL(GATK_HAPLOTYPE_CALLER.out, params.ref_fasta)
        SAMTOOLS_INDEX(LOFREQ_INDELQUAL.out,params.ref_fasta)
        LOFREQ_CALL(SAMTOOLS_INDEX.out)
        LOFREQ_FILTER(LOFREQ_CALL.out)

        // call_sv
        DELLY_CALL(GATK_MARK_DUPLICATES.out)
        BCFTOOLS_VIEW(DELLY_CALL.out)
        GATK_INDEX_FEATURE_FILE(BCFTOOLS_VIEW)
        //Enable this once a proper file with DR genes has been made:
        //$JAVA -Xmx64G -jar $GATK SelectVariants -V PREFIX.potentialSV.vcf.gz -O PREFIX.potentialSV.DRgenes.vcf.gz -L $RESOURCE_PATH/DRgenes.list

        // FIXME add the inputs
        // call_stats
        SAMTOOLS_STATS
        GATK_COLLECT_WGS_METRICS
        GATK_FLAG_STAT


        //TODO: Make sure that all of the stat files belong to the same sample (join operator)
        UTILS_SAMPLE_STATS

        UTILS_COHORT_STATS(UTILS_SAMPLE_STATS.out.collect())

}
