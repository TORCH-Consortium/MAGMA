FROM mambaorg/micromamba:1.4.9 AS base

#NOTE: The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER magma-tbprofiler-env.yml /tmp/magma-tbprofiler-env.yml

RUN micromamba install -y -f /tmp/magma-tbprofiler-env.yml -n base

RUN micromamba install -y -n base conda-forge::procps-ng  && micromamba clean -a -y

ENV PATH /opt/conda/bin:$PATH

FROM base AS tbdb


#NOTE: Just update the tb-profiler to the WHO-v2
#RUN tb-profiler update_tbdb --commit 22e9ceb6433d02071222af8dca506c6348604877
