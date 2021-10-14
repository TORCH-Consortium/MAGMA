nextflow.enable.dsl = 2


workflow SNP_ANALYSIS {

    // merge_select_snp
    GATK_SELECT_VARIANTS__SNP


    // merge_vqsr_snp
    GATK_VARIANT_RECALIBRATOR__SNP

    // merge_apply_vqsr_snp
    GATK_APPLY_VQSR__SNP
    GATK_SELECT_VARIANTS__SNP



}

workflow INDEL_ANALYSIS {
    // merge_select_indel
    GATK_SELECT_VARIANTS__INDEL

    // merge_vqsr_indel
    GATK_VARIANT_RECALIBRATOR__INDEL


    // merge_apply_vqsr_indel
    GATK_APPLY_VQSR__INDEL


}

workflow RESISTANCE_ANALYSIS {

    // merge_call_resistance
    TBPROFILER_VCF_PROFILE
    TBPROFILER_COLLATE

    // merge_call_resistance_lofreq
    TBPROFILER_VCF_PROFILE


}

workflow PHYLOGENY_ANALYSIS {
    // merge_phylogeny_prep_inccomplex
    GATK_SELECT_VARIANTS
    GATK_VARIANTS_TO_TABLE
    SNP_SITES
    SNP_DISTS

    // merge_iqtree_inccomplex
    IQTREE
}

workflow CLUSTER_PICKER_ANALYSIS {

    // merge_clusterpicker
    // 5/12 snp - including/excluding complex regions (4 trees)

}


workflow MERGE_WF {

    //TODO: collect joint stats - from utils

    //TODO: select samples based on that CSV

    //FIXME merge_combine
    GATK_COMBINE_GVCFS(params.vcf_name, FIXME_gvcfs, params.ref_fasta)


    // merge_genotype
    GATK_GENOTYPE_GVCFS(GATK_COMBINE_GVCFS.out, params.ref_fasta)

    // merge_snpeff_annotate
    GUNZIP
    SED
    SNPEFF
    SED
    //-----
    BGZIP
    //-----
    GATK_INDEX_FEATURE_FILE


    // merge_snp_indel_vcf
    GATK_MERGE_VCFS



    PHYLOGENY_ANALYSIS__INCCOMPLEX()
    PHYLOGENY_ANALYSIS__EXCOMPLEX()

    CLUSTER_PICKER_ANALYSIS__5SNP_INCCOMPLEX()
    CLUSTER_PICKER_ANALYSIS__12SNP_INCCOMPLEX()

    CLUSTER_PICKER_ANALYSIS__5SNP_EXCOMPLEX()
    CLUSTER_PICKER_ANALYSIS__12SNP_EXCOMPLEX()



}
