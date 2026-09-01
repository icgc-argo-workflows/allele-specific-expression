FROM continuumio/miniconda3:26.7.1-1

ARG PIPELINE_VERSION=dev

LABEL maintainer="adam.streck@iccb-cologne.org"
LABEL title="SchwarzLab ASE pipeline"
LABEL org.opencontainers.image.source="https://github.com/schwarzlab-ccb/pipeline-ase"
LABEL org.opencontainers.image.version="${PIPELINE_VERSION}"

COPY ase-environment.yaml conda-lock.yml /tmp/

RUN conda install --name base --channel conda-forge conda-lock=4.0.2 \
    && conda-lock install --prefix /opt/conda/envs/ase /tmp/conda-lock.yml \
    && conda clean --all --yes

ENV PATH=/opt/conda/envs/ase/bin:$PATH
