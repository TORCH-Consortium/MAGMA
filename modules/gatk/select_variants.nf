
$JAVA -Xmx64G -jar $GATK SelectVariants -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf.gz --select-type-to-include SNP -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz --remove-unused-alternates --exclude-non-variants


$JAVA -Xmx64G -jar $GATK SelectVariants -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf.gz --select-type-to-include INDEL -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_INDEL.vcf.gz --remove-unused-alternates --exclude-non-variants


$JAVA -Xmx64G -jar $GATK SelectVariants -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_inc-rRNA.vcf.gz -XL $RRNA -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_exc-rRNA.vcf.gz
