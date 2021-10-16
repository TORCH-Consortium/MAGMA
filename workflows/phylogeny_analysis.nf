
workflow PHYLOGENY_ANALYSIS {

    take:
    path(annotatedSnp)
    path(annotatedSnp)


    //----------
    // Including complex regions
    //----------

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


    // merge_phylogeny_prep_excomplex
    GATK_SELECT_VARIANTS('INDEL')
    GATK_VARIANTS_TO_TABLE
    SNP_SITES
    SNP_DISTS

    // merge_iqtree_excomplex
    IQTREE



}
