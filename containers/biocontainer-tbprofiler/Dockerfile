FROM quay.io/biocontainers/tb-profiler:6.0.0--pyhdfd78af_0 AS base

#FROM quay.io/biocontainers/tb-profiler:6.2.0--pyhdfd78af_0 AS base

FROM base AS whodb_v2

#NOTE: Just update the tb-profiler to the WHO-v2
RUN tb-profiler update_tbdb --commit 22e9ceb6433d02071222af8dca506c6348604877 --logging DEBUG
