FROM mambaorg/micromamba
MAINTAINER abhi18av@outlook.com

COPY --chown=$MAMBA_USER:$MAMBA_USER xbs-nf-env-2.yml /tmp/xbs-nf-env-2.yml

RUN micromamba install -y -f /tmp/xbs-nf-env-2.yml
