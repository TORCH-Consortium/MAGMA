tb-profiler vcf_profile --db $RESISTANCE_DB -d $OUT_DIR/resistance/$JOINT_NAME/XBS $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.filtered_SNP.RawIndels.vcf.gz


tb-profiler vcf_profile --lofreq_sample_name $i --db $RESISTANCE_DB -d $OUT_DIR/resistance/$JOINT_NAME/lofreq $OUT_DIR/minor_vars/$i.LoFreq.vcf.gz
