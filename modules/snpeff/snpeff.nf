
	 $JAVA -jar $SNPEFF -nostats -ud 40 Mycobacterium_tuberculosis_h37rv $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.vcf > $OUT_DIR/vcf/$JOINT_NAME/$JOINT_NAME.raw_variants.annotated.vcf



nextflow.enable.dsl = 2

params.results_dir = "${params.outdir}/snpeff"
params.save_mode = 'copy'
params.should_publish = true



process process_name {
    tag "something"
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    path(somefile)

    output:
    path("pattern"),  emit: "ch_output"

    script:

    """
    echo "nothing"
    """

    stub:

    """
    echo "nothing on stub"
    """

}