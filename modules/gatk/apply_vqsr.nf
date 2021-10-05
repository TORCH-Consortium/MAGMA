
         $JAVA -Xmx64G -jar $GATK ApplyVQSR -R $REFERENCE -V $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_SNP.vcf.gz -O $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP_inc-rRNA.vcf.gz --tranches-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.tranches --recal-file $OUT_DIR/vqsr/$JOINT_NAME/$JOINT_NAME.SNP.recal.vcf.gz --ts-filter-level 99.90 -AS --exclude-filtered -mode SNP
