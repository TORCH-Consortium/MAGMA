process UTILS_SAMPLE_QC {
    tag "joint_name: ${params.vcf_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(coverage), val(median_coverage_cutoff), path(wgsMetrics), path(flagStats), path(quanttbStats), path(ntmFraction)

    output:
    path("*.stats.tsv")
    path("QC_passed/*.g.vcf.gz")
    path("QC_passed/*.g.vcf.gz.tbi")


    shell:
//TODO: These could be added to directly the cohort CSV

    '''
$COVERAGE
$MEDIAN_COVERAGE_CUTOFF
$BREADTH_OF_COVERAGE
$BREADTH_OF_COVERAGE_CUTOFF
rel_abundance_threshold_met
ntm_fraction_threshold_met

    if [ "$COVERAGE" -ge "$MEDIAN_COVERAGE_CUTOFF" ] && [ 1 -eq "$(echo "$BREADTH_OF_COVERAGE >= $BREADTH_OF_COVERAGE_CUTOFF" | bc)" ] && rel_abundance_threshold_met $REL_ABUNDANCE_CUTOFF $REL_ABUNDANCE && ntm_fraction_threshold_met $NTM_FRACTION_CUTOFF $NTM_FRACTION;
    then
        echo "Adding ${i}"
        #echo "Adding ${i}, has a median coverage of ${COVERAGE} (${MEDIAN_COVERAGE_CUTOFF} required), breadth of coverage of ${BREADTH_OF_COVERAGE} (${BREADTH_OF_COVERAGE_CUTOFF} required), relative abundance of ${REL_ABUNDANCE} (${REL_ABUNDANCE_CUTOFF} required)"
        GVCFs+=" -V ${OUT_DIR}/gvcf/${i}.g.vcf.gz"
    else
        echo "Not adding ${i}, has a median coverage of ${COVERAGE} (${MEDIAN_COVERAGE_CUTOFF} minimum), breadth of coverage of ${BREADTH_OF_COVERAGE} (${BREADTH_OF_COVERAGE_CUTOFF} minimum), relative abundance of ${REL_ABUNDANCE} (${REL_ABUNDANCE_CUTOFF} minimum), NTM fraction of ${NTM_FRACTION} (${NTM_FRACTION_CUTOFF} maximum)"
    fi
    '''
