include { PREPARE_COHORT_VCF } from "./subworkflows/prepare_cohort_vcf.nf"
include { SNP_ANALYSIS } from "./subworkflows/snp_analysis.nf"
include { INDEL_ANALYSIS } from "./subworkflows/indel_analysis.nf"
include { GATK_MERGE_VCFS as GATK_MERGE_VCFS__INC } from "../modules/gatk/merge_vcfs.nf" addParams ( params.GATK_MERGE_VCFS )
include { GATK_MERGE_VCFS as GATK_MERGE_VCFS__EXC } from "../modules/gatk/merge_vcfs.nf"  addParams ( params.GATK_MERGE_VCFS )
include { MAJOR_VARIANT_ANALYSIS } from "./subworkflows/major_variant_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS__INCCOMPLEX } from "./subworkflows/phylogeny_analysis.nf"
include { PHYLOGENY_ANALYSIS as PHYLOGENY_ANALYSIS__EXCOMPLEX } from "./subworkflows/phylogeny_analysis.nf"
include { CLUSTER_ANALYSIS as CLUSTER_ANALYSIS__INCCOMPLEX } from "./subworkflows/cluster_analysis.nf"
include { CLUSTER_ANALYSIS as  CLUSTER_ANALYSIS__EXCOMPLEX } from "./subworkflows/cluster_analysis.nf"


workflow MERGE_WF {
    take:
        gvcf_ch
        reformatted_lofreq_vcfs_tuple_ch
        approved_samples_ch 

    main:

        //---------------------------------------------------------------------------------
        // Filter the approved samples
        //---------------------------------------------------------------------------------

        //NOTE: Reshape the flattened output of gvch_ch into the tuples of [sampleName, gvcf, gvcf.tbi]
        collated_gvcfs_ch = gvcf_ch
                                .flatten()
                                .collate(3)
                                .dump(tag:'MERGE_WF: collated_gvcfs_ch', pretty: true)
                                //.collectFile(name: "$params.outdir/collated_gvcfs_ch.txt")

        //NOTE: Join the fully approved samples with the gvcf channel 
        selected_gvcfs_ch = collated_gvcfs_ch.join(approved_samples_ch)
                                        .flatten()
                                        .dump(tag:'MERGE_WF: selected_gvcfs_ch', pretty: true)


        //NOTE: Filter only file type values and send to MERGE_WF
        //FIXME refactor the filtering logic NOT to rely upon the exact classnames
        filtered_selected_gvcfs_ch = selected_gvcfs_ch
                                        .filter { it.class != String }
                                        .collect()
                                        .dump(tag:'MERGE_WF: filtered_selected_gvcfs_ch', pretty: true)
                                        //.collectFile(name: "$params.outdir/selected_gvcfs_ch")


       //---------------------------------------------------------------------------------

        PREPARE_COHORT_VCF( filtered_selected_gvcfs_ch )

        SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        INDEL_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        merge_inc_vcf_ch = SNP_ANALYSIS.out.snp_inc_vcf_ch
                            .join(INDEL_ANALYSIS.out.indel_vcf_ch)
                            .dump(tag:'MERGE_WF: merge_inc_vcf_ch : ', pretty: true)

        // merge_snp_indel_vcf
        GATK_MERGE_VCFS__INC(merge_inc_vcf_ch)

        MAJOR_VARIANT_ANALYSIS(GATK_MERGE_VCFS__INC.out, reformatted_lofreq_vcfs_tuple_ch)


        //----------
        // Including complex regions
        //----------

        inccomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)])
                                                .ifEmpty([])
                                                .flatten()

        inccomplex_prefix_ch = Channel.of('ExDR.IncComplex')

        //NOTE: Both phylogenies should be excluding DR and excluding rRNA, then it is again filtered in two datasets one including complex regions and one excluding complex regions.
        //Ergo PHYLOGENY_...__INCCOMPLEX should take snp_exc_vcf_ch. Refer https://github.com/TORCH-Consortium/MAGMA/pull/114#discussion_r947732253
        PHYLOGENY_ANALYSIS__INCCOMPLEX(inccomplex_prefix_ch,
                                       inccomplex_exclude_interval_ref_ch,
                                       SNP_ANALYSIS.out.snp_exc_vcf_ch)

        CLUSTER_ANALYSIS__INCCOMPLEX(PHYLOGENY_ANALYSIS__INCCOMPLEX.out.snpsites_tree_tuple, inccomplex_prefix_ch)

        //----------
        // Excluding complex regions
        //----------

        excomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)],
                                                       [file(params.excluded_loci_list)])
                                                .ifEmpty([])
                                                .flatten()
                                                .dump(tag:'MERGE_WF: excomplex_exclude_interval_ref_ch', pretty: true)


        // excomplex_exclude_interval_ref_ch.view{ it -> "\n\n MAGMA-LOG MERGE_WF excomplex_exclude_interval_ref_ch: $it \n\n"}

        excomplex_prefix_ch = Channel.of('ExDR.ExComplex')


        PHYLOGENY_ANALYSIS__EXCOMPLEX(excomplex_prefix_ch,
                                       excomplex_exclude_interval_ref_ch,
                                       SNP_ANALYSIS.out.snp_exc_vcf_ch)

        CLUSTER_ANALYSIS__EXCOMPLEX(PHYLOGENY_ANALYSIS__EXCOMPLEX.out.snpsites_tree_tuple, excomplex_prefix_ch)


    emit:
        major_variants_results_ch =  MAJOR_VARIANT_ANALYSIS.out.major_variants_results_ch


}
