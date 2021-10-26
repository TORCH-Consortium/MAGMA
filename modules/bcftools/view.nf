process BCFTOOLS_VIEW {
    tag "${sampleName}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
    tuple val(sampleName), path(bcf)

    output:
    tuple val(sampleName), path(".*potentialSV.vcf.gz")

    shell:

    '''
    !{params.bcftools_path} view !{bcf} | grep '1/1' | sed 's/1\/1/1/g' | cut -f 1-3,5 -d "," > !{sampleName}.headerlessvcf.temp
    !{params.bcftools_path} view -h !{bcf} > !{sampleName}.header.temp

    cat !{sampleName}.header.temp !{sampleName}.headerlessvcf.temp \\
    | bgzip > !{sampleName}.potentialSV.vcf.gz

    '''

    stub:

    """
    touch ${sampleName}.potentialSV.vcf.gz
    """

}
