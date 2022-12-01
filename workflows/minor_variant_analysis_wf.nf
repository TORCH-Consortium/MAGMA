include { BCFTOOLS_MERGE } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE )
include { BGZIP } from "../modules/bgzip/bgzip.nf" addParams( params.BGZIP__MINOR_VARIANTS )
include { TBPROFILER_VCF_PROFILE__LOFREQ } from "../modules/tbprofiler/vcf_profile__lofreq.nf" addParams (params.TBPROFILER_VCF_PROFILE__LOFREQ)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__LOFREQ } from "../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__LOFREQ)
include { UTILS_MULTIPLE_INFECTION_FILTER } from "../modules/utils/multiple_infection_filter.nf" addParams (params.UTILS_MULTIPLE_INFECTION_FILTER)

workflow MINOR_VARIANT_ANALYSIS_WF {

    take:
        reformatted_lofreq_vcfs_tuple_ch

    main:

        vcfs_string_ch = reformatted_lofreq_vcfs_tuple_ch
                                .flatten()
                                .filter { it.extension  == "gz" }
                                .map { it -> it.name }
                                .reduce { a, b -> "$a $b " }
                                /* .view {"\n\n XBS-NF-LOG vcfs_filenames_ch : $it \n\n"} */

        BCFTOOLS_MERGE(vcfs_string_ch, reformatted_lofreq_vcfs_tuple_ch)

        // merge_call_resistance_lofreq
        BGZIP(BCFTOOLS_MERGE.out) 

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        //TBPROFILER minor variants
        TBPROFILER_VCF_PROFILE__LOFREQ(BGZIP.out, resistanceDb)

        //TBPROFILER major variants
        TBPROFILER_COLLATE__LOFREQ(params.vcf_name,
                                  TBPROFILER_VCF_PROFILE__LOFREQ.out.resistance_json.collect(),
                                  resistanceDb)

        UTILS_MULTIPLE_INFECTION_FILTER(TBPROFILER_COLLATE__LOFREQ.out.per_sample_results)

     emit: 
         approved_samples_ch = UTILS_MULTIPLE_INFECTION_FILTER.out.approved_samples_ch
         rejected_samples_ch = UTILS_MULTIPLE_INFECTION_FILTER.out.rejected_samples_ch

}
