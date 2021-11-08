include { GATK_COMBINE_GVCFS } from "../../modules/gatk/combine_gvcfs.nf" addParams ( params.GATK_COMBINE_GVCFS )
include { GATK_GENOTYPE_GVCFS } from "../../modules/gatk/genotyp_gvcfs.nf" addParams ( params.GATK_GENOTYPE_GVCFS )

workflow PREPARE_COHORT_VCF {
    take:
        cohort_gvcfs_ch


    main:
        // merge_combine
        GATK_COMBINE_GVCFS(params.vcf_name, cohort_gvcfs_ch, params.ref_fasta)


        // merge_genotype
        GATK_GENOTYPE_GVCFS(GATK_COMBINE_GVCFS.out, params.ref_fasta)

    /*
    // merge_snpeff_annotate
    SNPEFF(GATK_GENOTYPE_GVCFS)
    BGZIP(SNPEFF.out)
    GATK_INDEX_FEATURE_FILE(BGZIP.out)

    */
}
