include { GATK_COMBINE_GVCFS } from "../../modules/gatk/combine_gvcfs.nf" addParams ( params.GATK_COMBINE_GVCFS )
include { GATK_GENOTYPE_GVCFS } from "../../modules/gatk/genotype_gvcfs.nf" addParams ( params.GATK_GENOTYPE_GVCFS )
include { SNPEFF } from "../../modules/snpeff/snpeff.nf" addParams ( params.SNPEFF )
include { BGZIP } from "../../modules/bgzip/bgzip.nf" addParams( params.BGZIP )
include { GATK_INDEX_FEATURE_FILE__COHORT } from "../../modules/gatk/index_feature_file__cohort" addParams( params.GATK_INDEX_FEATURE_FILE__COHORT )


workflow PREPARE_COHORT_VCF {
    take:
        cohort_gvcfs_ch


    main:
        // merge_combine
        GATK_COMBINE_GVCFS(params.vcf_name, cohort_gvcfs_ch, params.ref_fasta)


        // merge_genotype
        GATK_GENOTYPE_GVCFS(GATK_COMBINE_GVCFS.out, params.ref_fasta)

        // merge_snpeff_annotate
        SNPEFF(GATK_GENOTYPE_GVCFS.out, params.ref_fasta)
        // BGZIP(SNPEFF.out)
        //TODO: Refactor to rely upon the singular INDEX_FILE_FEATURE module
        // GATK_INDEX_FEATURE_FILE__COHORT(BGZIP.out)

}
