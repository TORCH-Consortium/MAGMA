
$JAVA -Xmx64G -jar $GATK MergeVcfs -I $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_exc-rRNA.vcf.gz -I $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_INDEL.vcf.gz -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP.RawIndels.vcf.gz
