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

COPY main.py /tools/
ADD $fai_url /home/ubuntu
RUN awk -v OFS="\t" "{print \$1,\$2}" /home/ubuntu/GRCh38_Verily_v1.genome.fa.gz.fai | sort -V -s > /home/ubuntu/GRCh38_Verily_v1.genome.fa.gz.genome

# Add the bedgraph file from URL
# ADD $begraph_url /home/ubuntu

# Copy the bedgraph (local test only)
COPY k50.umap.bedgraph.gz /home/ubuntu/

RUN chmod 644 /home/ubuntu/*

WORKDIR /tools
USER ubuntu
ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
