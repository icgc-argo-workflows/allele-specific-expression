FROM broadinstitute/gatk:4.2.6.1

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/allele-specific-expression

RUN mkdir tools
WORKDIR /tools
ENV PATH="/tools:${PATH}:/tools/bcftools-1.9"

# Fix https://stackoverflow.com/questions/49582490/gpg-error-http-packages-cloud-google-com-apt-expkeysig-3746c208a7317b0f
RUN wget https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    && apt-key add apt-key.gpg
RUN apt-get update && apt-get install -y libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev

RUN cd /tools && \
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -vxjf bcftools-1.9.tar.bz2&& \
    cd bcftools-1.9 && \
    make

RUN groupadd -g 1000 ubuntu && \
   useradd -l -u 1000 -g ubuntu ubuntu && \
   install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

ARG fa_url=https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/rna-seq-references/GRCh38_Verily_v1.genome/GRCh38_Verily_v1.genome.fa.gz
ARG gzi_url=https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/rna-seq-references/GRCh38_Verily_v1.genome/GRCh38_Verily_v1.genome.fa.gz.gzi
ARG fai_url=https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/rna-seq-references/GRCh38_Verily_v1.genome/GRCh38_Verily_v1.genome.fa.gz.fai

ADD $fa_url /home/ubuntu
ADD $gzi_url /home/ubuntu
ADD $fai_url /home/ubuntu
RUN gatk CreateSequenceDictionary -R /home/ubuntu/GRCh38_Verily_v1.genome.fa.gz
RUN chmod 644 /home/ubuntu/*
# COPY GRCh38_Verily_v1.genome.fa.gz /home/ubuntu
# COPY GRCh38_Verily_v1.genome.fa.gz.fai /home/ubuntu
# COPY GRCh38_Verily_v1.genome.fa.gz.gzi /home/ubuntu


USER ubuntu
ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
