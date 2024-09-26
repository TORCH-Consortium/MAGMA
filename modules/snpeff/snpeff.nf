process SNPEFF {
    tag "${joint_name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(joint_name), path(rawJointVariantsFile)
        path(ref_fasta)

    output:
        tuple val(joint_name), path("*raw_variants.annotated.vcf")


    script:

        """
        rename_vcf_chrom.py --vcf ${rawJointVariantsFile}  --source ${params.ref_fasta_basename} --target 'Chromosome' --outfile temp.chromosomes.vcf

        ${params.snpeff_path} -nostats ${params.arguments} temp.chromosomes.vcf > temp.annotated.vcf

        rename_vcf_chrom.py --target ${params.ref_fasta_basename} --source 'Chromosome' --vcf temp.annotated.vcf --outfile ${joint_name}.temp.vcf

        cp ${joint_name}.temp.vcf  ${joint_name}.raw_variants.annotated.vcf
        """

    stub:

        """
        touch ${joint_name}.raw_variants.annotated.vcf
        """

}
