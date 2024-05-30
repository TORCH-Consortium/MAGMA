FROM quay.io/biocontainers/tb-profiler:6.2.1--pyhdfd78af_0 AS base

FROM base AS tbprofiler

#NOTE: Just update the tb-profiler databaes to rely upon the relevant branch.

# WHO-v2 specific tag https://github.com/jodyphelan/tbdb/releases/tag/who-v2-strict
# COMMIT bdace1f82d948ce0001e1dade6eb93d2da9c47e5 

# MAGMA branch 
#RUN tb-profiler update_tbdb --branch MAGMA --logging DEBUG

RUN tb-profiler update_tbdb --commit 30f8bc37df15affa378ebbfbd3e1eb4c5903056e --logging DEBUG
