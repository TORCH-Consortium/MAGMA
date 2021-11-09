include { GATK_SELECT_VARIANTS as  GATK_SELECT_VARIANTS__SNP } from "../../modules/gatk/select_variants.nf" addParams ( params.GATK_SELECT_VARIANTS__SNP )
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR__SNP } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__SNP )
include { GATK_APPLY_VQSR as GATK_APPLY_VQSR__SNP } from "../../modules/gatk/apply_vqsr.nf" addParams ( params.GATK_APPLY_VQSR__SNP )

workflow SNP_ANALYSIS {

    take:
        cohort_vcf_and_index_ch

    main:

        // merge_select_snp
        GATK_SELECT_VARIANTS__SNP('SNP', cohort_vcf_and_index_ch, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict] )


        // merge_vqsr_snp

        arg_files_ch = Channel.of(["coll2014,known=false,training=true,truth=true,prior=15.0", file(params.coll2014_vcf), file(params.coll2014_vcf_tbi)],
                              ["coll2018,known=false,training=true,truth=true,prior=15.0", file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)],
                              ["Napier2020,known=false,training=true,truth=true,prior=15.0", file(params.napier2020_vcf), file(params.napier2020_vcf_tbi)],
                              ["Benavente2015,known=true,training=false,truth=false,prior=5.0", file(params.benavente2015_vcf), file(params.benavente2015_vcf_tbi)])
            .ifEmpty([])
            .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1], it[2] ] : [] }
            .flatten()
            // .view()


        args_ch = arg_files_ch
            .filter { it.class == org.codehaus.groovy.runtime.GStringImpl }
            .reduce { a, b -> "$a --resource:$b " }
            .ifEmpty("")
            // .view()


        resources_files_ch = arg_files_ch
            .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
            .filter {  it.getExtension()  == "gz" }
            .collect()
            .ifEmpty([])
            // .view()

        resources_file_indexes_ch = arg_files_ch
            .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
            .filter {  it.getExtension()  == "tbi" }
            .collect()
            .ifEmpty([])
            // .view()



        GATK_VARIANT_RECALIBRATOR__SNP('SNP',
                                    GATK_SELECT_VARIANTS__SNP.out,
                                    args_ch,
                                    resources_files_ch,
                                    resources_file_indexes_ch,
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )


        vqsr_ch = GATK_SELECT_VARIANTS__SNP.out.variantsVcfTuple
            .join(GATK_VARIANT_RECALIBRATOR__SNP.out.recalVcfTuple)
            .join(GATK_VARIANT_RECALIBRATOR__SNP.out.tranchesFile)

        // merge_apply_vqsr_snp
        GATK_APPLY_VQSR__SNP('SNP',
                            vqsr_ch,
                            params.ref_fasta,
                            [params.ref_fasta_fai, params.ref_fasta_dict])

/*
        GATK_SELECT_VARIANTS__EXCLUSION__SNP('SNP', GATK_APPLY_VQSR__SNP.out, params.rrna_file)


    emit:
        GATK_SELECT_VARIANTS__EXCLUSION__SNP.out

*/
}
