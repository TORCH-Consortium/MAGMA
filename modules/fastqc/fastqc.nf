process FASTQC {
    tag "${sampleName}"
    label 'cpu_2_memory_2'
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        tuple val(sampleName), val(bamRgString), path(sampleReads)

    output:
        path('*fastqc*')


    script:

        """
        ${params.fastqc_path} \\
            ${sampleReads} \\
            -t ${task.cpus}
        """

    stub:
        """
        echo "${params.fastqc_path} \\
                ${sampleReads} \\
                -t ${task.cpus}"

        touch ${sampleName}.html

        touch ${sampleName}.zip
        """
}
