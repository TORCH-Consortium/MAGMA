
workflow INDEL_ANALYSIS {

    // INDEL analysis is experimental XBS_merge#L164

    // merge_select_indel
    GATK_SELECT_VARIANTS('INDEL', GATK_INDEX_FEATURE_FILE.out )

    // merge_vqsr_indel

    arg_files_ch = Channel.of(["walker2015,known=true,training=true,truth=true,prior=15.0", file("Walker2015.UVPapproved.rRNAexcluded.vcf.gz")],
                              //FIXME XBS_merge#L169
                              //["zeng2018,known=true,training=true,truth=true,prior=15.0", file("Walker2015.UVPapproved.rRNAexcluded.vcf.gz") ]
                              ["coll2018,known=false,training=true,truth=true,prior=15.0", file("Coll2018.UVPapproved.rRNAexcluded.vcf.gz")])
    .ifEmpty([])
    .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1] ] : [] }
    .flatten()


args_ch = arg_files_ch
    .filter { it.class == org.codehaus.groovy.runtime.GStringImpl }
    .reduce { a, b -> "$a --resource:$b " }
    .ifEmpty("")
    .view()


files_ch = arg_files_ch
    .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
    .collect()
    .ifEmpty([])
    .view()




    GATK_VARIANT_RECALIBRATOR__INDEL('INDEL', GATK_SELECT_VARIANTS__INDEL.out, args_ch, files_ch )


    // merge_apply_vqsr_snp
    GATK_APPLY_VQSR__INDEL('INDEL', GATK_SELECT_VARIANTS__INDEL.out, GATK_VARIANT_RECALIBRATOR__INDEL.out.recalVcf, params.ref_fasta)


    emit:
    //NOTE: This is supposed to be temporary output XBS_merge#L189
    GATK_SELECT_VARIANTS.out

}
