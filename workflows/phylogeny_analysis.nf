
workflow PHYLOGENY_ANALYSIS {

    take:
    path(annotatedSnp)
    path(annotatedSnp)


    //----------
    // Including complex regions
    //----------


    arg_files_ch = Channel.of([file("Coll2018.vcf.gz")])
        .ifEmpty([])
        .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1] ] : [] }
        .flatten()


    args_ch = arg_files_ch
        .filter { it.class == org.codehaus.groovy.runtime.GStringImpl }
        .reduce { a, b -> "$a -XL:$b " }
        .ifEmpty("")
        .view()


    files_ch = arg_files_ch
        .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
        .collect()
        .ifEmpty([])
        .view()


    // merge_phylogeny_prep_inccomplex
    GATK_SELECT_VARIANTS('SNP')
    GATK_VARIANTS_TO_TABLE
    SNP_SITES
    SNP_DISTS

    // merge_iqtree_inccomplex
    IQTREE


    //----------
    // Excluding complex regions
    //----------

    arg_files_ch = Channel.of([file("UVP_List_of_Excluded_loci.list")
                               file("Coll2018.vcf.gz")])
        .ifEmpty([])
        .map { it -> it != [] ? [ "${it[0]} ${it[1].getName()}", it[1] ] : [] }
        .flatten()


    args_ch = arg_files_ch
        .filter { it.class == org.codehaus.groovy.runtime.GStringImpl }
        .reduce { a, b -> "$a -XL:$b " }
        .ifEmpty("")
        .view()


    files_ch = arg_files_ch
        .filter { it.class != org.codehaus.groovy.runtime.GStringImpl }
        .collect()
        .ifEmpty([])
        .view()



    // merge_phylogeny_prep_excomplex
    GATK_SELECT_VARIANTS('INDEL')
    GATK_VARIANTS_TO_TABLE
    SNP_SITES
    SNP_DISTS

    // merge_iqtree_excomplex
    IQTREE



}
