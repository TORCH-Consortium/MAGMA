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
include { UTILS_MERGE_COHORT_STATS } from "../modules/utils/merge_cohort_stats.nf" addParams ( params.UTILS_MERGE_COHORT_STATS )


workflow MERGE_WF {
    take:
        gvcf_ch
        reformatted_lofreq_vcfs_tuple_ch
        cohort_stats_tsv
        approved_samples_ch
        rejected_samples_ch

    main:

        //---------------------------------------------------------------------------------
        // Filter the approved samples
        //---------------------------------------------------------------------------------

        //NOTE: Read the approved_samples tsv file and isolate the names of the approved samples
        approved_samples_minor_variants_ch = approved_samples_ch
                                .splitCsv(header: false, skip: 1, sep: '\t' )
                                .map { row -> [ row.first() ] }
                                .collect()
                                .dump(tag:'MERGE_WF: approved_samples_minor_variants_ch', pretty: true)
                                /* .view {"\n\n XBS-NF-LOG approved_samples_minor_variants_ch : $it \n\n"} */

        //NOTE: Reshape the flattened output of gvch_ch into the tuples of [sampleName, gvcf, gvcf.tbi]
        collated_gvcfs_ch = gvcf_ch
                                .flatten()
                                .collate(3)
                                .dump(tag:'MERGE_WF: collated_gvcfs_ch', pretty: true)
                                /* .view {"\n\n XBS-NF-LOG collated_gvcfs_ch : $it \n\n"} */
                                //.collectFile(name: "$params.outdir/collated_gvcfs_ch.txt")



        //FIXME: Refactor this to emit two different files and use only the approved samples
        //NOTE: Use the stats file for the entire cohort (from CALL_WF)
        // and filter out the samples which pass all thresholds
        approved_call_wf_samples_ch = cohort_stats_tsv
                                .splitCsv(header: false, skip: 1, sep: '\t' )
                                .map { row -> [
                                        row.first(),           // SAMPLE
                                        row.last().toInteger() // ALL_THRESHOLDS_MET
                                        ]
                                    }
                                .filter { it[1] == 1} // Filter out samples which meet all the thresholds
                                .map { [ it[0] ] }
                                .dump(tag:'MERGE_WF: approved_call_wf_samples_ch', pretty: true)

        /* approved_call_wf_samples_ch */
        /*         .collect() */
        /*         .dump(tag:'approved_call_wf_samples_ch.collect()') */
        /*         .view {"\n\n XBS-NF-LOG approved_call_wf_samples_ch.collect() : $it \n\n"} */ 

        //NOTE: Join the approved samples from MINOR_VARIANT_ANALYSIS_WF and CALL_WF
        fully_approved_samples_ch = approved_samples_minor_variants_ch
                                        .join(approved_call_wf_samples_ch)
                                        .flatten()
                                        .dump(tag:'MERGE_WF: fully_approved_samples_ch', pretty: true)
                                        /* .view {"\n\n XBS-NF-LOG fully_approved_samples_ch : $it \n\n"} */
                                        //.collect()
                                        //.collectFile(name: "$params.outdir/approved_samples_ch.txt") 


        //NOTE: Join the fully approved samples with the gvcf channel 
        selected_gvcfs_ch = collated_gvcfs_ch.join(fully_approved_samples_ch)
                                        .flatten()
                                        .dump(tag:'MERGE_WF: selected_gvcfs_ch', pretty: true)

        //NOTE: Filter only file type values and send to MERGE_WF
        filtered_selected_gvcfs_ch = selected_gvcfs_ch
                                        .filter { it -> { 
                                                            (it.class.name  == "sun.nio.fs.UnixPath") 
                                                            || (it.class.name == "nextflow.cloud.azure.nio.AzPath") 
                                                            || (it.class.name == "com.upplication.s3fs.S3Path") 
                                                            || (it.class.name == "com.google.cloud.storage.contrib.nio.CloudStoragePath") 
                                                    } }
                                        .collect()
                                        .dump(tag:'MERGE_WF: filtered_selected_gvcfs_ch', pretty: true)
                                        //.collectFile(name: "$params.outdir/selected_gvcfs_ch")


       //---------------------------------------------------------------------------------

        UTILS_MERGE_COHORT_STATS()

       //---------------------------------------------------------------------------------

        PREPARE_COHORT_VCF(selected_gvcfs_ch)

        SNP_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        INDEL_ANALYSIS(PREPARE_COHORT_VCF.out.cohort_vcf_and_index_ch)

        merge_inc_vcf_ch = SNP_ANALYSIS.out.snp_inc_vcf_ch
                            .join(INDEL_ANALYSIS.out.indel_vcf_ch)
                            .dump(tag:'MERGE_WF: merge_inc_vcf_ch : ', pretty: true)

        // merge_snp_indel_vcf
        GATK_MERGE_VCFS__INC(merge_inc_vcf_ch)

        MAJOR_VARIANT_ANALYSIS(GATK_MERGE_VCFS__INC.out, reformatted_lofreq_vcf_ch)


        //----------
        // Including complex regions
        //----------

        inccomplex_exclude_interval_ref_ch = Channel.of([file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)])
                                                .ifEmpty([])
                                                .flatten()

        inccomplex_prefix_ch = Channel.of('ExDR.IncComplex')

        //NOTE: Both phylogenies should be excluding DR and excluding rRNA, then it is again filtered in two datasets one including complex regions and one excluding complex regions.
        //Ergo PHYLOGENY_...__INCCOMPLEX should take snp_exc_vcf_ch. Refer https://github.com/TORCH-Consortium/xbs-nf/pull/114#discussion_r947732253
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


        // excomplex_exclude_interval_ref_ch.view{ it -> "\n\n XBS-NF-LOG MERGE_WF excomplex_exclude_interval_ref_ch: $it \n\n"}

        excomplex_prefix_ch = Channel.of('ExDR.ExComplex')


        PHYLOGENY_ANALYSIS__EXCOMPLEX(excomplex_prefix_ch,
                                       excomplex_exclude_interval_ref_ch,
                                       SNP_ANALYSIS.out.snp_exc_vcf_ch)

        CLUSTER_ANALYSIS__EXCOMPLEX(PHYLOGENY_ANALYSIS__EXCOMPLEX.out.snpsites_tree_tuple, excomplex_prefix_ch)


    emit:
        major_variants_results_ch =  MAJOR_VARIANT_ANALYSIS.out.major_variants_results_ch


}
