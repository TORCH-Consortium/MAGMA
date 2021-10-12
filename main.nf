nextflow.enable.dsl = 2


//================================================================================
// Derive file names and location from the params.yaml
//================================================================================

//================================================================================
// Include sub-workflows and (soft) override workflow-level parameters
//================================================================================


workflow MAP_WF {

    // - The RG is derived from CSV fields
    // - Accomodates both single/paired ends
    BWA_MEM
}


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


workflow MERGE_WF {

    // collect joint stats - from utils

    // merge_combine
    GATK_COMBINE_GVCFS


    // merge_genotype
    GATK_GENOTYPE_GVCFS

    // merge_snpeff_annotate
    GUNZIP
    SED
    SNPEFF
    SED
    GATK_INDEX_FEATURE_FILE


    // merge_select_snp
    GATK_SELECT_VARIANTS

    // merge_select_indel
    GATK_SELECT_VARIANTS


    // merge_vqsr_snp
    GATK_VARIANT_RECALIBRATOR

    // merge_apply_vqsr_snp
    GATK_APPLY_VQSR
    GATK_SELECT_VARIANTS

    // merge_vqsr_indel
    GATK_VARIANT_RECALIBRATOR


    // merge_apply_vqsr_indel
    GATK_APPLY_VQSR


    // merge_snp_indel_vcf
    GATK_MERGE_VCFS


    // merge_call_resistance
    TBPROFILER_VCF_PROFILE
    TBPROFILER_COLLATE

    // merge_call_resistance_lofreq
    TBPROFILER_VCF_PROFILE


    // merge_phylogeny_prep_inccomplex
    GATK_SELECT_VARIANTS
    GATK_VARIANTS_TO_TABLE
    SNP_SITES
    SNP_DISTS

    // merge_iqtree_inccomplex
    IQTREE

    // merge_phylogeny_prep_excomplex
    GATK_SELECT_VARIANTS
    GATK_VARIANTS_TO_TABLE
    SNP_SITES
    SNP_DISTS



    // merge_iqtree_excomplex
    IQTREE

    // merge_clusterpicker
    // 5/12 snp - including/excluding complex regions (4 trees)



}

//================================================================================
// Prepare channels
//================================================================================

//================================================================================
// Main workflow
//================================================================================

workflow {

    FASTQC()
    MAP_WF()

    QUANTTB(MAP_WF.out)
    CALL_WF(MAP_WF.out)

    MERGE_WF(QUANTTB.out, CALL_WF.out)

}
