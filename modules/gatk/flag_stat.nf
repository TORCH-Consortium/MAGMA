
$JAVA -Xmx64G -jar $GATK FlagStat -R $REFERENCE -I $OUT_DIR/mapped/$SAMPLE_ID.recal_reads.bam > $OUT_DIR/stats/$SAMPLE_ID.recal_reads.FlagStat
