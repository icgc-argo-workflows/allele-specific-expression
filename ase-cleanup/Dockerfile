FROM python 

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/allele-specific-expression

RUN apt-get update && apt-get install -y bedtools
RUN pip install pandas scipy statsmodels

ENV PATH="/tools:${PATH}"

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

ARG fai_url=https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/rna-seq-references/GRCh38_Verily_v1.genome/GRCh38_Verily_v1.genome.fa.gz.fai
ARG begraph_url=https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bedgraph.gz 

COPY *.py /tools/
WORKDIR /tools
ADD $fai_url /tools/
RUN awk -v OFS="\t" "{print \$1,\$2}" GRCh38_Verily_v1.genome.fa.gz.fai | sort -s > GRCh38_Verily_v1.genome.fa.gz.genome
ADD $begraph_url /tools/ 

# local test only
# COPY k50.umap.bedgraph.gz /tools/k50.umap.bedgraph.gz

USER ubuntu
ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
