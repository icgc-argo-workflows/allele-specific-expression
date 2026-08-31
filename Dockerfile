FROM continuumio/miniconda3

ARG PIPELINE_VERSION=dev

LABEL maintainer="adam.streck@gmail.com"
LABEL title="SchwarzLab ASE pipeline"
LABEL org.opencontainers.image.source="https://github.com/schwarzlab-ccb/pipeline-ase"
LABEL org.opencontainers.image.version="${PIPELINE_VERSION}"

COPY ase-environment.yaml /tmp/ase-environment.yaml

RUN conda env create --file /tmp/ase-environment.yaml \
    && conda clean --all --yes

ENV PATH=/opt/conda/envs/ase/bin:$PATH
