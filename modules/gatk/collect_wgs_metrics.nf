
$JAVA -Xmx64G -jar $GATK CollectWgsMetrics -R $REFERENCE -I $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam -O $OUT_DIR/stats/$SAMPLE_ID.recal_reads.WgsMetrics --READ_LENGTH 0 --COVERAGE_CAP 10000 --COUNT_UNPAIRED
