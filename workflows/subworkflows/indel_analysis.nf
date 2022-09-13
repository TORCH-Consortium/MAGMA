include { GATK_SELECT_VARIANTS as  GATK_SELECT_VARIANTS__INDEL } from "../../modules/gatk/select_variants.nf" addParams ( params.GATK_SELECT_VARIANTS__INDEL  )
include { GATK_VARIANT_RECALIBRATOR as GATK_VARIANT_RECALIBRATOR__INDEL  } from "../../modules/gatk/variant_recalibrator.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__INDEL  )
include { GATK_APPLY_VQSR as GATK_APPLY_VQSR__INDEL  } from "../../modules/gatk/apply_vqsr.nf" addParams ( params.GATK_APPLY_VQSR__INDEL  )
include { GATK_SELECT_VARIANTS__EXCLUSION as  GATK_SELECT_VARIANTS__EXCLUSION__INDEL  } from "../../modules/gatk/select_variants__exclusion.nf" addParams ( params.GATK_SELECT_VARIANTS__EXCLUSION__INDEL  )
include { OPTIMIZE_VARIANT_RECALIBRATION } from "./optimize_variant_recalibration.nf" addParams ( params.GATK_VARIANT_RECALIBRATOR__INDEL )



//NOTE: INDEL analysis is experimental XBS_merge#L164
workflow INDEL_ANALYSIS {

    take:
        cohort_vcf_and_index_ch

    main:

        // merge_select_indel
        GATK_SELECT_VARIANTS__INDEL('INDEL',
                                    'raw',
                                    cohort_vcf_and_index_ch,
                                    "",
                                    [],
                                    [],
                                    params.ref_fasta,
                                    [params.ref_fasta_fai, params.ref_fasta_dict] )

    /*
        //NOTE: This section isn't used as of now, XBS_merge#L189
        // merge_vqsr_indel

        arg_files_ch = Channel.of(["walker2015,known=true,training=true,truth=true,prior=15.0", file(params.walker2015_vcf), file(params.walker2015_vcf_tbi)],
                                ["walker2015,known=true,training=true,truth=true,prior=15.0", file(params.zeng2018_vcf), file(params.zeng2018_vcf_tbi)],
                                ["coll2018,known=false,training=true,truth=true,prior=15.0", file(params.coll2018_vcf), file(params.coll2018_vcf_tbi)])
            .ifEmpty([])
            .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1], it[2] ] : [] }
            .flatten()


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


        //FIXME XBS_merge#L169
        if(params.optimize_variant_recalibration) {

            OPTIMIZE_VARIANT_RECALIBRATION('INDEL',
                                        GATK_SELECT_VARIANTS__INDEL.out.variantsVcfTuple,
                                        args_ch,
                                        resources_files_ch,
                                        resources_file_indexes_ch)

            vqsr_ch = OPTIMIZE_VARIANT_RECALIBRATION.out.optimized_vqsr_ch

        } else {

            GATK_VARIANT_RECALIBRATOR__INDEL('INDEL',
                                    " -an DP -an AS_MQ ",
                                        GATK_SELECT_VARIANTS__INDEL.out.variantsVcfTuple,
                                        args_ch,
                                        resources_files_ch,
                                        resources_file_indexes_ch,
                                        params.ref_fasta,
                                        [params.ref_fasta_fai, params.ref_fasta_dict] )

            vqsr_ch = GATK_SELECT_VARIANTS__INDEL.out.variantsVcfTuple
                .join(GATK_VARIANT_RECALIBRATOR__INDEL.out.recalVcfTuple)
                .join(GATK_VARIANT_RECALIBRATOR__INDEL.out.tranchesFile)

        }


        // merge_apply_vqsr_snp
        GATK_APPLY_VQSR__INDEL('INDEL',
                            vqsr_ch,
                            params.ref_fasta,
                            [params.ref_fasta_fai, params.ref_fasta_dict])

    */

    emit:
        //NOTE: This is supposed to be temporary output XBS_merge#L189 and should be replaced with the output of GATK_APPLY_VQSR__INDEL
        indel_vcf_ch = GATK_SELECT_VARIANTS__INDEL.out.variantsVcfTuple

        // indel_vcf_ch = GATK_APPLY_VQSR__INDEL.out.filteredVcfTuple

}
