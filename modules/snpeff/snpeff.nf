process SNPEFF {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(joint_name), path(rawJointVariantsFile)
    path(ref_fasta)

    output:
    tuple val(joint_name), path(".*annotated.vcf")


    script:

    """
    gunzip -c ${rawJointVariantsFile} > ${joint_name}.raw_variants.vcf

    sed -i 's/^${ref_fasta.getBaseName()}/Chromosome/g' ${joint_name}.raw_variants.vcf

    ${params.snpeff_path} \\
        ${params.arguments} \\
        Mycobacterium_tuberculosis_h37rv \\
        ${rawJointVariantsFile} \\
    > ${joint_name}.raw_variants.annotated.vcf

    rm ${joint_name}.raw_variants.vcf

    sed -i 's/^Chromosome/${ref_fasta.getBaseName()}/g' ${joint_name}.raw_variants.annotated.vcf

    """

    stub:

    """
    touch ${joint_name}.raw_variants.annotated.vcf
    """

}
