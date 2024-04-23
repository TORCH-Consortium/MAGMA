FROM quay.io/biocontainers/tb-profiler:6.2.0--pyhdfd78af_0 AS base

FROM base AS tbprofiler

#NOTE: Just update the tb-profiler to the WHO-v2 https://github.com/jodyphelan/tbdb/releases/tag/who-v2-strict
RUN tb-profiler update_tbdb --commit bdace1f82d948ce0001e1dade6eb93d2da9c47e5 --logging DEBUG
