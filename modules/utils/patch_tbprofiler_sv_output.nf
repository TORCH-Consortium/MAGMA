process UTILS_PATCH_TBPROFILER_SV_OUTPUT {
    publishDir params.results_dir, mode: params.save_mode, enabled: params.should_publish

    input:
        path vcf_files
        tuple val(joint_name), path(joint_vcf_index), path(joint_vcf)

    shell:

        '''

        tbdb_path=$(tb-profiler list_db | grep -o '/[^ ]*')

        patch_tbprofiler_sv_output.py \\
            --delly_vcf \\
            --ismapper_vcf \\
            --bed_file "$tbdb_path/tbdb.bed" \\
            --existing_json_file \\
            --cleaned_json_file \\
            --output_vcf

        '''


}
