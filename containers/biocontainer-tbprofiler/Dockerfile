FROM quay.io/biocontainers/tb-profiler:6.3.0--pyhdfd78af_0 AS base

FROM base AS tbprofiler

#NOTE: Just update the tb-profiler databaes to rely upon the relevant branch.

# WHO-v2 specific tag https://github.com/jodyphelan/tbdb/releases/tag/who-v2-strict
# COMMIT bdace1f82d948ce0001e1dade6eb93d2da9c47e5 

# MAGMA branch 
#RUN tb-profiler update_tbdb --branch MAGMA --logging DEBUG

RUN tb-profiler update_tbdb --commit 70c0f203d9ca69aaf8efb39c8d6865be3306ee92 --logging DEBUG
