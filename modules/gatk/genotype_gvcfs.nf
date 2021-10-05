
$JAVA -Xmx64G -jar $GATK GenotypeGVCFs -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.combined.vcf.gz -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation --sample-ploidy 1
