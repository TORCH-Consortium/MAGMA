process TBPROFILER_COLLATE {
    publishdir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    //FIXME This is only relevant for the
    beforeScript "source ${params.conda_dir}/etc/profile.d/conda.sh"
    conda "${params.tb_profiler_venv}"
    afterScript "conda deactivate"


    input:
    tuple val(joint_name), path(resistanceDb)

    output:
    path("*.XBS.resistance*")

    script:
    """
    tb-profiler collate --db ${resistanceDb} -p ${joint_name}.XBS.resistance.txt
    """

    stub:
    """
    touch ${joint_name}.XBS.resistance.txt
    """
}


