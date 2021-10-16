nextflow.enable.dsl = 2


workflow MERGE_WF {

    //TODO: collect joint stats - from utils

    //TODO: select samples based on that CSV

    PREPARE_COHORT_VCF_WF()

    SNP_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)
    INDEL_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)

    // merge_snp_indel_vcf
    GATK_MERGE_VCFS(SNP_ANALYSIS_WF.out.filteredVCF, INDEL_ANALYSIS_WF.out.filteredVCF)

    RESISTANCE_ANALYSIS()


    PHYLOGENY_ANALYSIS__INCCOMPLEX()
    PHYLOGENY_ANALYSIS__EXCOMPLEX()

    CLUSTER_ANALYSIS(PHYLOGENY_ANALYSIS__INCCOMPLEX.out, PHYLOGENY_ANALYSIS__EXCOMPLEX.out)

}
