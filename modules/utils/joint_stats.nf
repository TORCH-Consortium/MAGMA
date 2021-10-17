

echo -e "SAMPLE\tAVG_INSERT_SIZE\tMAPPED_%\tRAW_TOTAL_SEQS\tAVERAGE_QUALITY\tQUANTTB_RELATIVE_ABUNDANCE\tQUANTTB_DEPTH\tMEAN_COVERAGE\tSD_COVERAGE\tMEDIAN_COVERAGE\tMAD_COVERAGE\tPCT_EXC_ADAPTER\tPCT_EXC_MAPQ\tPCT_EXC_DUPE\tPCT_EXC_UNPAIRED\tPCT_EXC_BASEQ\tPCT_EXC_OVERLAP\tPCT_EXC_CAPPED\tPCT_EXC_TOTAL\tPCT_1X\tPCT_5X\tPCT_10X\tPCT_30X\tPCT_50X\tPCT_100X\tNTM_FRACTION" > ${OUT_DIR}/stats/$JOINT_NAME.summary.tsv

read -ra SAMPLE_LIST <<< "$SAMPLES"
for i in "${SAMPLE_LIST[@]}"
do
    COVERAGE=$(cat ${OUT_DIR}/stats/${i}.recal_reads.WgsMetrics | grep "^4411532" | cut -f 4)
    BREADTH_OF_COVERAGE=$(cat ${OUT_DIR}/stats/${i}.recal_reads.WgsMetrics | grep "^4411532" | cut -f 14)
    REL_ABUNDANCE=$(tail -n +2 ${OUT_DIR}/quanttb/${i}.qtb.txt | cut -d , -f 4 | paste -s -d ';')
    DEPTH=$(tail -n +2 ${OUT_DIR}/quanttb/${i}.qtb.txt | cut -d , -f 5 | paste -s -d ';')
    TOTAL_SEQS=$(cat ${OUT_DIR}/stats/${i}.recal_reads.SamtoolStats | grep "insert size average" | cut -f 3)
    MAPPED_P=$(cat ${OUT_DIR}/stats/${i}.recal_reads.FlagStat | grep "mapped (" | cut -f 2 -d "(" | cut -f 1 -d "%")
    INS_SIZE=$(cat ${OUT_DIR}/stats/${i}.recal_reads.SamtoolStats | grep "raw total sequences" | cut -f 3)
    AVG_QUAL=$(cat ${OUT_DIR}/stats/${i}.recal_reads.SamtoolStats | grep "average quality" | cut -f 3)
    WGS_METR=$(cat ${OUT_DIR}/stats/${i}.recal_reads.WgsMetrics | grep "^4411532" | cut -f 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20,22,27)
    NTM_FRACTION=$(cat ${OUT_DIR}/stats/${i}.potential_NTM_fraction.txt)

    echo -e "${i}\t${TOTAL_SEQS}\t${MAPPED_P}\t${INS_SIZE}\t${AVG_QUAL}\t${REL_ABUNDANCE}\t${DEPTH}\t${WGS_METR}\t${NTM_FRACTION}" >> ${OUT_DIR}/stats/$JOINT_NAME.summary.tsv

    if [ "$COVERAGE" -ge "$MEDIAN_COVERAGE_CUTOFF" ] && [ 1 -eq "$(echo "$BREADTH_OF_COVERAGE >= $BREADTH_OF_COVERAGE_CUTOFF" | bc)" ] && rel_abundance_threshold_met $REL_ABUNDANCE_CUTOFF $REL_ABUNDANCE && ntm_fraction_threshold_met $NTM_FRACTION_CUTOFF $NTM_FRACTION;
    then
        echo "Adding ${i}"
        #echo "Adding ${i}, has a median coverage of ${COVERAGE} (${MEDIAN_COVERAGE_CUTOFF} required), breadth of coverage of ${BREADTH_OF_COVERAGE} (${BREADTH_OF_COVERAGE_CUTOFF} required), relative abundance of ${REL_ABUNDANCE} (${REL_ABUNDANCE_CUTOFF} required)"
        GVCFs+=" -V ${OUT_DIR}/gvcf/${i}.g.vcf.gz"
    else
        echo "Not adding ${i}, has a median coverage of ${COVERAGE} (${MEDIAN_COVERAGE_CUTOFF} minimum), breadth of coverage of ${BREADTH_OF_COVERAGE} (${BREADTH_OF_COVERAGE_CUTOFF} minimum), relative abundance of ${REL_ABUNDANCE} (${REL_ABUNDANCE_CUTOFF} minimum), NTM fraction of ${NTM_FRACTION} (${NTM_FRACTION_CUTOFF} maximum)"
    fi;
done
