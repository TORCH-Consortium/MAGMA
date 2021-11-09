include { PREPARE_COHORT_VCF } from "./subworkflows/prepare_cohort_vcf.nf"
include { SNP_ANALYSIS } from "./subworkflows/snp_analysis.nf"
// include { INDEL_ANALYSIS } from "./subworkflows/indel_analysis.nf"

workflow MERGE_WF {
    take:
        selected_gvcfs_ch


    main:
        PREPARE_COHORT_VCF(selected_gvcfs_ch)

        SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

    /*
        INDEL_ANALYSIS_WF(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        // merge_snp_indel_vcf
        GATK_MERGE_VCFS(SNP_ANALYSIS_WF.out.filteredVCF, INDEL_ANALYSIS_WF.out.filteredVCF)

        RESISTANCE_ANALYSIS()


        PHYLOGENY_ANALYSIS__INCCOMPLEX()
        PHYLOGENY_ANALYSIS__EXCOMPLEX()

        CLUSTER_ANALYSIS(PHYLOGENY_ANALYSIS__INCCOMPLEX.out, PHYLOGENY_ANALYSIS__EXCOMPLEX.out)
    */

}
