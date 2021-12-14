process TBPROFILER_LOAD_LIBRARY {

    input:
        path(resistanceDb)

    output:
        path(resistanceDb)

    script:

        if (!workflow.container) {

            """
            echo "Load the library here as no BioContainer based container is used"

            cd ${resistanceDb}

            ${params.tbprofiler_path} load_library ${resistanceDb.name}
            """

        } else {

            """
            echo "Do NOT load the library here as the no BioContainer based container is used"
            """
        }

    stub:
        """
        echo "This is a container based run => ${workflow.container}"
        echo "${params.tbprofiler_path} load_library ${resistanceDb.name}"
        """

}
