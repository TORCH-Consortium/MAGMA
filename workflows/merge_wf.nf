include { PREPARE_COHORT_VCF } from "./subworkflows/prepare_cohort_vcf.nf"
include { SNP_ANALYSIS } from "./subworkflows/snp_analysis.nf"
include { INDEL_ANALYSIS } from "./subworkflows/indel_analysis.nf"
include { GATK_MERGE_VCFS } from "../modules/gatk/merge_vcfs.nf" addParams ( params.GATK_MERGE_VCFS )
include { RESISTANCE_ANALYSIS } from "./subworkflows/resistance_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS_INCCOMPLEX } from "./subworkflows/phylogeny_analysis.nf"
// include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS_EXCOMPLEX } from "./subworkflows/phylogeny_analysis.nf"


workflow MERGE_WF {
    take:
        selected_gvcfs_ch
        lofreq_vcf_ch


    main:
        PREPARE_COHORT_VCF(selected_gvcfs_ch)

        SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        INDEL_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        merge_vcf_ch = (SNP_ANALYSIS.out.snp_vcf_ch).join(INDEL_ANALYSIS.out.indel_vcf_ch)

        // merge_snp_indel_vcf
        GATK_MERGE_VCFS(merge_vcf_ch)

        RESISTANCE_ANALYSIS(GATK_MERGE_VCFS.out, lofreq_vcf_ch)

        PHYLOGENY_ANALYSIS__INCCOMPLEX('ExDR.IncComplex',
                                    SNP_ANALYSIS.out.snp_vcf_ch,
                                    params.rrna)

    /*
        PHYLOGENY_ANALYSIS__EXCOMPLEX()

        CLUSTER_ANALYSIS(PHYLOGENY_ANALYSIS__INCCOMPLEX.out, PHYLOGENY_ANALYSIS__EXCOMPLEX.out)
    */

}
