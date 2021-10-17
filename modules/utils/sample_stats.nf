
// SAMPLE    AVG_INSERT_SIZE      MAPPED_%    RAW_TOTAL_SEQS    AVERAGE_QUALITY    QUANTTB_RELATIVE_ABUNDANCE     RELATIVE_ABUNDANCE_THRESHOLD_MET    QUANTTB_DEPTH    MEAN_COVERAGE    SD_COVERAGE       NTM_FRACTION_THRESHOLD_MET MEDIAN_COVERAGE    MAD_COVERAGE    PCT_EXC_ADAPTER    PCT_EXC_MAPQ    PCT_EXC_DUPE    PCT_EXC_UNPAIRED    PCT_EXC_BASEQ    PCT_EXC_OVERLAP    PCT_EXC_CAPPED    PCT_EXC_TOTAL    PCT_1X    PCT_5X    PCT_10X    PCT_30X    PCT_50X    PCT_100X    NTM_FRACTION
// ${i}      ${TOTAL_SEQS}      ${MAPPED_P}    ${INS_SIZE}    ${AVG_QUAL}          ${REL_ABUNDANCE}                 $(rel_abundance_threshold_met)       ${DEPTH}       ${WGS_METR}      ${NTM_FRACTION}   $(ntm_fraction_threshold_met)


process UTILS_SAMPLE_STATS {

    input:
    tuple val(sampleName), path(samtoolsStats), path(wgsMetrics), path(flagStats), path(quanttbStats), path(ntmFraction)

    output:
    path("*.stats.tsv")


    shell:

    '''
// this function checks whther the relative abundance threshold is met for one of the strains
function rel_abundance_threshold_met () {
    IFS=';'
    local i
    for i in ${@:2}
    do
        if [ 1 -eq "$(echo "$i >= $1" | bc)" ]
        then
            IFS=' '
            echo 0
            return 0
        fi
    done
    IFS=' '
    echo 1
    return 1
}


// this function checks whether the relative abundance threshold is met for one of the strains
function ntm_fraction_threshold_met () {
    IFS=';'
    local i
    for i in ${@:2}
    do
        if [ 1 -eq "$(echo "$i <= $1" | bc)" ]
        then
            IFS=' '
            echo 0
            return 0
        fi
    done
    IFS=' '
    echo 1
    return 1
}



    COVERAGE=$(cat !{wgsMetrics} | grep "^4411532" | cut -f 4)
    BREADTH_OF_COVERAGE=$(cat !{wgsMetrics} | grep "^4411532" | cut -f 14)
    REL_ABUNDANCE=$(tail -n +2 !{quanttbStats} | cut -d , -f 4 | paste -s -d ';')
    DEPTH=$(tail -n +2 !{quanttbStats} | cut -d , -f 5 | paste -s -d ';')
    TOTAL_SEQS=$(cat !{samtoolsStats} | grep "insert size average" | cut -f 3)
    MAPPED_P=$(cat !{flagStats} | grep "mapped (" | cut -f 2 -d "(" | cut -f 1 -d "%")
    INS_SIZE=$(cat !{samtoolsStats} | grep "raw total sequences" | cut -f 3)
    AVG_QUAL=$(cat !{samtoolsStats} | grep "average quality" | cut -f 3)
    WGS_METR=$(cat !{wgsMetrics} | grep "^4411532" | cut -f 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,22,27)
    NTM_FRACTION=$(cat !{ntmFraction})

    echo -e "${i}\t${TOTAL_SEQS}\t${MAPPED_P}\t${INS_SIZE}\t${AVG_QUAL}\t${REL_ABUNDANCE}\t$(rel_abundance_threshold_met $REL_ABUNDANCE_CUTOFF $REL_ABUNDANCE)\t${DEPTH}\t${WGS_METR}\t${NTM_FRACTION}\t$(ntm_fraction_threshold_met $NTM_FRACTION_CUTOFF $NTM_FRACTION)" > !{sampleName}.stats.tsv
    '''
}
