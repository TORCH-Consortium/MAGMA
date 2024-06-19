include { GATK_COMBINE_GVCFS } from "../../modules/gatk/combine_gvcfs.nf" addParams ( params.GATK_COMBINE_GVCFS )
include { GATK_GENOTYPE_GVCFS } from "../../modules/gatk/genotype_gvcfs.nf" addParams ( params.GATK_GENOTYPE_GVCFS )
include { SNPEFF } from "../../modules/snpeff/snpeff.nf" addParams ( params.SNPEFF )
include { BGZIP } from "../../modules/bgzip/bgzip.nf" addParams( params.BGZIP )
include { GATK_INDEX_FEATURE_FILE as GATK_INDEX_FEATURE_FILE__COHORT } from "../../modules/gatk/index_feature_file" addParams( params.GATK_INDEX_FEATURE_FILE__COHORT )


workflow PREPARE_COHORT_VCF {
    take:
        cohort_gvcfs_ch


    main:
        //FIXME save the string to an intermediate file
        gvcfs_string_ch = cohort_gvcfs_ch
            .flatten()
            .filter {  file(it).getExtension()  == "gz" }
            .map { it -> file(it).name }
            .reduce { a, b -> "$a --variant $b " }
            .dump(tag: "PREPARE_COHORT_VCF gvcfs_string_ch:", pretty: true)


        if (params.use_ref_gvcf) {
            refExitRifGvcf = file(params.ref_gvcf, checkIfExists: true)
        } else {
            refExitRifGvcf = []
        }

        if (params.use_ref_gvcf) {
            refExitRifGvcfTbi = file(params.ref_gvcf_tbi, checkIfExists: true)
        } else {
            refExitRifGvcfTbi = []
        }

        // merge_combine
        GATK_COMBINE_GVCFS(params.vcf_name,
                        gvcfs_string_ch,
                        cohort_gvcfs_ch,
                        params.ref_fasta,
                        refExitRifGvcf,
                        refExitRifGvcfTbi,
                        [params.ref_fasta_fai, params.ref_fasta_dict])


        // merge_genotype
        GATK_GENOTYPE_GVCFS(GATK_COMBINE_GVCFS.out, params.ref_fasta, [params.ref_fasta_fai, params.ref_fasta_dict])

        // merge_snpeff_annotate
        SNPEFF(GATK_GENOTYPE_GVCFS.out, params.ref_fasta)
        BGZIP(SNPEFF.out)
        GATK_INDEX_FEATURE_FILE__COHORT(BGZIP.out, '')

    emit:
        cohort_vcf_and_index_ch = GATK_INDEX_FEATURE_FILE__COHORT.out.sample_vcf_tuple
}
