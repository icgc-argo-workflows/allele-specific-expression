FROM python 

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/allele-specific-expression

RUN apt-get update && apt-get install -y bedtools
RUN pip install pandas scipy statsmodels

ENV PATH="/tools:${PATH}"

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

COPY *.py /tools/
RUN cd /tools
RUN wget https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bedgraph.gz
RUN wget https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa.gz.fai
RUN awk -v OFS=\"\\t\" \"{print \\\$1,\\\$2}\" GRCh38_hla_decoy_ebv.fa.gz.fai | sort -s > GRCh38_hla_decoy_ebv.fa.gz.genome
WORKDIR /tools    


USER ubuntu
ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
