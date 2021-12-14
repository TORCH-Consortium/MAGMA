process TBPROFILER_LOAD_LIBRARY {

    input:
        path(resistanceDb)

    script:
        """
        cd ${resistanceDb}

        ${params.tbprofiler_path} load_library ${resistanceDb.name}
        """

    stub:
        """
        echo "${params.tbprofiler_path} load_library ${resistanceDb.name}"
        """

}
