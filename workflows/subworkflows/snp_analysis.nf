include { GATK_SELECT_VARIANTS as  GATK_SELECT_VARIANTS__SNP } from "../../modules/gatk/select_variants.nf" addParams ( params.GATK_SELECT_VARIANTS )

workflow SNP_ANALYSIS {

    take:
        cohort_vcf_and_index_ch

    main:


        gvcfs_ch = cohort_vcf_and_index_ch.view()


        // merge_select_snp
        // GATK_SELECT_VARIANTS__SNP('SNP', cohort_vcf_and_index_ch, [params.ref_fasta, params.ref_fasta_fai] )


        // merge_vqsr_snp

    // arg_files_ch = Channel.of([

    //     "coll2014,known=false,training=true,truth=true,prior=15.0", file("${projectDir}/resources/truth/Coll2014.UVPapproved.rRNAexcluded.vcf.gz")],


    //                               // ["coll2014,known=false,training=true,truth=true,prior=15.0", file(params.coll2014_vcf)],
    //                               // ["coll2018,known=false,training=true,truth=true,prior=15.0", file(params.coll2018_vcf)],
    //                               // ["Napier2020,known=false,training=true,truth=true,prior=15.0", file(params.napier2020_vcf)],
    //                               // ["Benavente2015,known=true,training=false,truth=false,prior=5.0", file(params.benavente2015_vcf)]
    // )
    //         .ifEmpty([])
    //         .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1] ] : [] }
    //         .flatten()


        // args_ch = arg_files_ch
        //     .filter { it.class == org.codehaus.groovy.runtime.GStringImpl }
        //     .reduce { a, b -> "$a --resource:$b " }
        //     .ifEmpty("")
        //     .view()


        // files_ch = arg_files_ch
        //     .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
        //     .collect()
        //     .ifEmpty([])
        //     .view()


/*
        GATK_VARIANT_RECALIBRATOR__SNP('SNP', GATK_SELECT_VARIANTS__SNP.out, args_ch, files_ch )



        // merge_apply_vqsr_snp
        GATK_APPLY_VQSR__SNP('SNP', GATK_SELECT_VARIANTS__SNP.out, GATK_VARIANT_RECALIBRATOR__SNP.out.recalVcf, params.ref_fasta)

        GATK_SELECT_VARIANTS__EXCLUSION__SNP('SNP', GATK_APPLY_VQSR__SNP.out, params.rrna_file)


    emit:
        GATK_SELECT_VARIANTS__EXCLUSION__SNP.out

*/
}
