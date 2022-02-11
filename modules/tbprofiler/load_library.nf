process TBPROFILER_LOAD_LIBRARY {

    input:
        path(resistanceDb)

    output:
        path(resistanceDb)

    script:

        if (!workflow.container && resistanceDb) {

            """
            echo "A container ${workflow.container} is NOT used: TRUE"
            echo "resistanceDB status: ${resistanceDb}"

            echo "Load the library"

            cd ${resistanceDb}

            ${params.tbprofiler_path} load_library ${resistanceDb.name}
            """

        } else {

            """
            echo "A container ${workflow.container} is used: TRUE"
            echo "resistanceDB status: ${resistanceDb.name}"


            echo "Do NOT load the library"
            """
        }

    stub:
        """
        echo "This is a container based run => ${workflow.container}"
        echo "${params.tbprofiler_path} load_library ${resistanceDb.name}"
        """

}
