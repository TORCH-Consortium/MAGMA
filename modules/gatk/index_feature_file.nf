	 gunzip -c $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf.gz > $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf
	 sed -i 's/^NC-000962-3-H37Rv/Chromosome/g' $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf
	 $JAVA -jar $SNPEFF -nostats -ud 40 Mycobacterium_tuberculosis_h37rv $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf > $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf
	 rm $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf
	 sed -i 's/^Chromosome/NC-000962-3-H37Rv/g' $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf
	 $BGZIP -f $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf
     $JAVA -Xmx64G -jar $GATK IndexFeatureFile -I $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf.gz
