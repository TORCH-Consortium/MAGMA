process TBPROFILER_COLLATE {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    //FIXME This is only relevant for the PBS configuration and should be mentioned there after Nextflow v21.10 release
    // beforeScript "source ${params.conda_dir}/etc/profile.d/conda.sh"
    // conda "${params.tb_profiler_venv}"
    // afterScript "conda deactivate"


    input:
    tuple val(joint_name), path("results/*")
    path(resistanceDb)

    output:
    path("*.XBS.resistance*")

    script:
    def optionalDb  = resistanceDb ? "--db ${resistanceDb}" : ""

    """
    ${params.tbprofiler_path} collate \\
        ${optionalDb} \\
        -p ${joint_name}.XBS.resistance
    """

    stub:
    """
    touch ${joint_name}.XBS.resistance.txt
    """
}


