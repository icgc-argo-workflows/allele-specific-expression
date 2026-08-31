FROM continuumio/miniconda3

LABEL maintainer="adam.streck@gmail.com"
LABEL title="SchwarzLab ASE pipeline"

COPY ase-environment.yaml /tmp/ase-environment.yaml

RUN conda env create --file /tmp/ase-environment.yaml \
    && conda clean --all --yes

ENV PATH=/opt/conda/envs/ase/bin:$PATH
