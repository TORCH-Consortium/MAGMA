include { PREPARE_COHORT_VCF } from "../modules/samtools/merge.nf" addParams ( params.SAMTOOLS_MERGE )

workflow MERGE_WF {
    take:
        path(cohort_stats_file)



    //TODO: select samples based on that CSV - publish files for both (passing/failing) samples
    // PREPARE_COHORT_VCF(cohort_stats_file)

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
