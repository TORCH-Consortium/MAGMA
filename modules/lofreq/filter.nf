process LOFREQ_FILTER {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(vcf)
    path(ref_fasta)

    output:
    tuple val(sampleName), path("*.Filtered_AF.vcf")

    script:

    """
	lofreq filter \\
        ${arguments} \\
        -i $OUT_DIR/minor_vars/${vcf}  \\
    > ${sampleName}.Filtered_AF.vcf
    """

    stub:

    """
    echo "lofreq filter \\
        ${arguments} \\
        -i $OUT_DIR/minor_vars/${vcf}  \\
        ${sampleName}.Filtered_AF.vcf"

    touch ${sampleName}.Filtered_AF.vcf

    """

}
