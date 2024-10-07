process UTILS_PATCH_TBPROFILER_SV_OUTPUT {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path("delly_and_ismapper_concat_vcfs/*") //path vcf_files
        path("tbprofiler_json") //tuple val(joint_name), path(joint_vcf_index), path(joint_vcf)

    shell:

        '''
        tbdb_path=$(tb-profiler list_db | grep -o '/[^ ]*')


        patch_tbprofiler_sv_output.py \\
            --bed_file "$tbdb_path/tbdb.bed" \\
            --output_vcf \\
            --existing_json_file \\
            --cleaned_json_file

        '''

}
