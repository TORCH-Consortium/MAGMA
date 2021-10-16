nextflow.enable.dsl = 2


workflow MERGE_WF {

    //TODO: collect joint stats - from utils

    //TODO: select samples based on that CSV

    PREPARE_COHORT_VCF_WF()
    SNP_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)
    INDEL_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)

    // merge_snp_indel_vcf
    GATK_MERGE_VCFS


    PHYLOGENY_ANALYSIS
    PHYLOGENY_ANALYSIS

    CLUSTER_ANALYSIS()

}
