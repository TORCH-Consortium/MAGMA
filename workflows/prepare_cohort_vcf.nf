workflow PREPARE_COHORT_VCF {

    //FIXME merge_combine
    GATK_COMBINE_GVCFS(params.vcf_name, FIXME_gvcfs, params.ref_fasta)


    // merge_genotype
    GATK_GENOTYPE_GVCFS(GATK_COMBINE_GVCFS.out, params.ref_fasta)

    // merge_snpeff_annotate
    SNPEFF(GATK_GENOTYPE_GVCFS)
    BGZIP(SNPEFF.out)
    GATK_INDEX_FEATURE_FILE(BGZIP.out)


}
