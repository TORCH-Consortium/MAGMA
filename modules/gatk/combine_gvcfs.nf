
	check_exit $JAVA -Xmx64G -jar $GATK CombineGVCFs -R $REFERENCE -G StandardAnnotation -G AS_StandardAnnotation $GVCFs -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.combined.vcf.gz
