
	check_exit $SAMTOOLS stats -F DUP,SUPPLEMENTARY,SECONDARY,UNMAP,QCFAIL $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam -r $REFERENCE > $OUT_DIR/stats/$SAMPLE_ID.recal_reads.SamtoolStats
