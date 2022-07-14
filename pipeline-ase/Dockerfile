FROM continuumio/miniconda3

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/allele-specific-expression

RUN conda config --add channels default
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda


RUN pip install --no-cache-dir --upgrade pip
RUN pip install pandas numpy scipy statsmodels pyensembl

RUN conda install -c bioconda gatk4 -y
RUN conda install -c bioconda samtools=1.14 -y
RUN conda install -c bioconda bedtools -y
RUN conda install -c bioconda bcftools -y
RUN conda install -c bioconda pyvcf -y
RUN conda install -c bioconda matplotlib -y

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
