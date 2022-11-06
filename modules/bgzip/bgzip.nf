process BGZIP {
    tag "${name}"
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(name), path(annotatedVcf)

    output:
        tuple val(name), path("*.gz")

    script:

        """
        ${params.bgzip_path} ${params.arguments} ${annotatedVcf}
        """

    stub:

        """
        touch ${name}.annotated.vcf.gz
        """

}
