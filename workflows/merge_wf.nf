include { PREPARE_COHORT_VCF } from "./subworkflows/prepare_cohort_vcf.nf"

workflow MERGE_WF {
    take:
        selected_gvcfs_ch


    main:
        PREPARE_COHORT_VCF(selected_gvcfs_ch)

    /*
    SNP_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)
    INDEL_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)

    // merge_snp_indel_vcf
    GATK_MERGE_VCFS(SNP_ANALYSIS_WF.out.filteredVCF, INDEL_ANALYSIS_WF.out.filteredVCF)

    RESISTANCE_ANALYSIS()


    PHYLOGENY_ANALYSIS__INCCOMPLEX()
    PHYLOGENY_ANALYSIS__EXCOMPLEX()

    CLUSTER_ANALYSIS(PHYLOGENY_ANALYSIS__INCCOMPLEX.out, PHYLOGENY_ANALYSIS__EXCOMPLEX.out)
    */

}
