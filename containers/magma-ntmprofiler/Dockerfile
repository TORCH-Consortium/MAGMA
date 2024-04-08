FROM mambaorg/micromamba:1.5.8 AS base

#NOTE: The conda env file has been copied via the build script
COPY --chown=$MAMBA_USER:$MAMBA_USER magma-ntmprofiler-env.yml /tmp/magma-ntmprofiler-env.yml

RUN micromamba install -y -f /tmp/magma-ntmprofiler-env.yml -n base

#RUN micromamba install -y -n base conda-forge::procps-ng  && micromamba clean -a -y

#FROM base AS ntmprofiler

#ENV PATH /opt/conda/bin:$PATH

#RUN ntm-profiler update_db --commit 93c979b4c3825c591e5bceb55f86538e10d556d2

