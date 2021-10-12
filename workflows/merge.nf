nextflow.enable.dsl = 2

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
