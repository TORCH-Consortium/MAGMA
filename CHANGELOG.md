# CHANGELOG FOR THE MAGMA PIPELINE VERSIONS
<!-- https://keepachangelog.com/en/1.1.0/ -->


## v2.0.0

1. Update `tb-profiler` to `v6.2.1` along with the `WHO mutations catalog v2` with the database version `30f8bc37df15affa378ebbfbd3e1eb4c5903056e`
 
2. Addition of `ntm-profiler` with the database FIXME

3. Two new `conda-envs` and containers have been added to the pipeline, relying directly on upstream biocontainers for both `ntm-profiler` and `tb-profiler`

4. Changes to the directory structure

5. Addition of mixed-infections resistance summary script

6. Addition of `low_memory` profile to accommodate lower end infrastructure

7. Updates to the names of parameters for triggering partial/optional workflows


## v1.1.1

Created a parallel workflow for mapping without using the strict seed lenght for use in the structural variant workflow.

Updated TBProfiler to version 5.0.0 and recreated the resistance database to work with the the new version

Updated the summarize resistance script to include the structural variants in the excel output


## v1.0.0

Initial release of the pipeline
