#FROM mambaorg/micromamba:1.4.9 AS base

#NOTE: The conda env file has been copied via the build script
#COPY --chown=$MAMBA_USER:$MAMBA_USER magma-ntmprofiler-env.yml /tmp/magma-ntmprofiler-env.yml

#RUN micromamba install -y -f /tmp/magma-ntmprofiler-env.yml -n base

#RUN micromamba install -y -n base conda-forge::procps-ng  && micromamba clean -a -y

#FROM base AS ntmdb

#ENV PATH /opt/conda/bin:$PATH

#RUN ntm-profiler update_db --commit 93c979b4c3825c591e5bceb55f86538e10d556d2 --logging DEBUG

#=================================================
# FAILED ATTEMPT TO PATCH EXISTING BIOCONTAINER 
# https://app.clickup.com/t/86941gv3t?comment=90120030159229
# https://app.clickup.com/t/86941gv3t?comment=90120033586266
#=================================================

#NOTE: Make sure the version of ntm-profiler version is correct
FROM quay.io/biocontainers/ntm-profiler:0.4.0--pyhdfd78af_0

#NOTE: Just update the ntm-profiler to the latest version of the database
RUN ntm-profiler update_db  --logging DEBUG

