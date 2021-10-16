nextflow.enable.dsl = 2


workflow MERGE_WF {

    //TODO: collect joint stats - from utils

    //TODO: select samples based on that CSV


    if [ "$COVERAGE" -ge "$MEDIAN_COVERAGE_CUTOFF" ] && [ 1 -eq "$(echo "$BREADTH_OF_COVERAGE >= $BREADTH_OF_COVERAGE_CUTOFF" | bc)" ] && rel_abundance_threshold_met $REL_ABUNDANCE_CUTOFF $REL_ABUNDANCE && ntm_fraction_threshold_met $NTM_FRACTION_CUTOFF $NTM_FRACTION;
    then
        echo "Adding ${i}"
        #echo "Adding ${i}, has a median coverage of ${COVERAGE} (${MEDIAN_COVERAGE_CUTOFF} required), breadth of coverage of ${BREADTH_OF_COVERAGE} (${BREADTH_OF_COVERAGE_CUTOFF} required), relative abundance of ${REL_ABUNDANCE} (${REL_ABUNDANCE_CUTOFF} required)"
        GVCFs+=" -V ${OUT_DIR}/gvcf/${i}.g.vcf.gz"
    else
        echo "Not adding ${i}, has a median coverage of ${COVERAGE} (${MEDIAN_COVERAGE_CUTOFF} minimum), breadth of coverage of ${BREADTH_OF_COVERAGE} (${BREADTH_OF_COVERAGE_CUTOFF} minimum), relative abundance of ${REL_ABUNDANCE} (${REL_ABUNDANCE_CUTOFF} minimum), NTM fraction of ${NTM_FRACTION} (${NTM_FRACTION_CUTOFF} maximum)"
    fi;



    PREPARE_COHORT_VCF_WF()

    SNP_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)
    INDEL_ANALYSIS_WF(PREPARE_COHORT_VCF_WF.out)

    // merge_snp_indel_vcf
    GATK_MERGE_VCFS(SNP_ANALYSIS_WF.out.filteredVCF, INDEL_ANALYSIS_WF.out.filteredVCF)

    RESISTANCE_ANALYSIS()


    PHYLOGENY_ANALYSIS__INCCOMPLEX()
    PHYLOGENY_ANALYSIS__EXCOMPLEX()

    CLUSTER_ANALYSIS(PHYLOGENY_ANALYSIS__INCCOMPLEX.out, PHYLOGENY_ANALYSIS__EXCOMPLEX.out)

}
