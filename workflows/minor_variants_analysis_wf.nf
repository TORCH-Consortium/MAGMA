include { BCFTOOLS_MERGE as BCFTOOLS_MERGE__LOFREQ } from "../modules/bcftools/merge.nf" addParams ( params.BCFTOOLS_MERGE__LOFREQ )
include { BGZIP } from "../modules/bgzip/bgzip.nf" addParams( params.BGZIP__MINOR_VARIANTS )
include { TBPROFILER_VCF_PROFILE__LOFREQ } from "../modules/tbprofiler/vcf_profile__lofreq.nf" addParams (params.TBPROFILER_VCF_PROFILE__LOFREQ)
include { TBPROFILER_COLLATE as TBPROFILER_COLLATE__LOFREQ } from "../modules/tbprofiler/collate.nf" addParams (params.TBPROFILER_COLLATE__LOFREQ)
include { UTILS_MULTIPLE_INFECTION_FILTER } from "../modules/utils/multiple_infection_filter.nf" addParams (params.UTILS_MULTIPLE_INFECTION_FILTER)
include { UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION } from "../modules/utils/summarize_resistance_results_mixed_infection.nf" addParams (params.UTILS_SUMMARIZE_RESISTANCE_RESULTS_MIXED_INFECTION)

workflow MINOR_VARIANTS_ANALYSIS_WF {

    take:
        reformatted_lofreq_vcfs_tuple_ch

    main:

//FIXME save the string to an intermediate file
        vcfs_string_ch = reformatted_lofreq_vcfs_tuple_ch
                                .flatten()
                                .filter { it.extension  == "gz" }
                                .map { it -> it.name }
                                //.reduce { a, b -> "$a $b " }
                                //.dump(tag:'MINOR_VARIANT_WF: vcfs_string_ch', pretty: true)

        // merge_call_resistance_lofreq
        //NOTE: Samples implicitly get filtered here if they don't have any identified variants
	      vcfs_file = vcfs_string_ch.collectFile(name: 'minor_variant_vcfs.txt', newLine: true)

        BCFTOOLS_MERGE__LOFREQ(vcfs_file, reformatted_lofreq_vcfs_tuple_ch)

        def resistanceDb =  params.resistance_db != "NONE" ?  params.resistance_db : []

        TBPROFILER_VCF_PROFILE__LOFREQ(BCFTOOLS_MERGE__LOFREQ.out, resistanceDb)

        TBPROFILER_COLLATE__LOFREQ(params.vcf_name,
                                  TBPROFILER_VCF_PROFILE__LOFREQ.out.resistance_json.collect(),
                                  resistanceDb)

        UTILS_MULTIPLE_INFECTION_FILTER(TBPROFILER_COLLATE__LOFREQ.out.per_sample_results)


     emit: 
         approved_samples_ch = UTILS_MULTIPLE_INFECTION_FILTER.out.approved_samples
         rejected_samples_ch = UTILS_MULTIPLE_INFECTION_FILTER.out.rejected_samples
         minor_variants_results_ch = TBPROFILER_COLLATE__LOFREQ.out.per_sample_results

}
