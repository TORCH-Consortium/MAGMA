include { PREPARE_COHORT_VCF } from "./subworkflows/prepare_cohort_vcf.nf"
include { SNP_ANALYSIS } from "./subworkflows/snp_analysis.nf"
include { INDEL_ANALYSIS } from "./subworkflows/indel_analysis.nf"
include { GATK_MERGE_VCFS } from "../modules/gatk/merge_vcfs.nf" addParams ( params.GATK_MERGE_VCFS )
include { RESISTANCE_ANALYSIS } from "./subworkflows/resistance_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS__INCCOMPLEX } from "./subworkflows/phylogeny_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS__EXCOMPLEX } from "./subworkflows/phylogeny_analysis.nf"
include { CLUSTER_ANALYSIS as CLUSTER_ANALYSIS__INCCOMPLEX } from "./subworkflows/cluster_analysis.nf"
include { CLUSTER_ANALYSIS as  CLUSTER_ANALYSIS__EXCOMPLEX } from "./subworkflows/cluster_analysis.nf"


workflow MERGE_WF {
    take:
        selected_gvcfs_ch
        lofreq_vcf_ch


    main:
        PREPARE_COHORT_VCF(selected_gvcfs_ch)

        SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        INDEL_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        merge_vcf_ch = (SNP_ANALYSIS.out.snp_vcf_ch).join(INDEL_ANALYSIS.out.indel_vcf_ch)

        merge_vcf_ch.view( it -> "\n\n XBS-NF-LOG MERGE_WF merge_vcf_ch: $merge_vcf_ch \n\n")

        // merge_snp_indel_vcf
        GATK_MERGE_VCFS(merge_vcf_ch)

        RESISTANCE_ANALYSIS(GATK_MERGE_VCFS.out, lofreq_vcf_ch)


        //----------
        // Including complex regions
        //----------

        inccomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)])
                                                .ifEmpty([])
                                                .flatten()

        inccomplex_prefix_ch = Channel.of('ExDR.IncComplex')


        PHYLOGENY_ANALYSIS__INCCOMPLEX(inccomplex_prefix_ch,
                                       inccomplex_exclude_interval_ref_ch,
                                       SNP_ANALYSIS.out.snp_vcf_ch)

        CLUSTER_ANALYSIS__INCCOMPLEX(PHYLOGENY_ANALYSIS__INCCOMPLEX.out.snpsites_tree_tuple, inccomplex_prefix_ch)

        //----------
        // Excluding complex regions
        //----------

        excomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)],
                                                       [file(params.excluded_loci_list)])
                                                .ifEmpty([])
                                                .flatten()


        excomplex_exclude_interval_ref_ch.view( it -> "\n\n XBS-NF-LOG MERGE_WF excomplex_exclude_interval_ref_ch: $excomplex_exclude_interval_ref_ch \n\n")

        excomplex_prefix_ch = Channel.of('ExDR.ExComplex')


        PHYLOGENY_ANALYSIS__EXCOMPLEX(excomplex_prefix_ch,
                                       excomplex_exclude_interval_ref_ch,
                                       SNP_ANALYSIS.out.snp_vcf_ch)

        CLUSTER_ANALYSIS__EXCOMPLEX(PHYLOGENY_ANALYSIS__EXCOMPLEX.out.snpsites_tree_tuple, excomplex_prefix_ch)


}
