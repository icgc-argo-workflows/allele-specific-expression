FROM python

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/allele-specific-expression

# setuptools 57 due to https://stackoverflow.com/questions/69100275/error-while-downloading-the-requirements-using-pip-install-setup-command-use-2
RUN pip install setuptools==57 
RUN pip install pyvcf pandas scipy statsmodels matplotlib pyensembl

ENV PATH="/tools:${PATH}"

RUN groupadd -g 1000 ubuntu && \
    useradd -l -u 1000 -g ubuntu ubuntu && \
    install -d -m 0755 -o ubuntu -g ubuntu /home/ubuntu

ARG genome_url=https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/rna-seq-references/GRCh38_Verily_v1.annotation/gencode.v40.chr_patch_hapl_scaff.annotation.gtf

COPY *.py /tools/

# Download the genome annotation
ADD $genome_url /home/ubuntu
RUN /tools/createGtfDB.py GRCh38 /home/ubuntu/gencode.v40.chr_patch_hapl_scaff.annotation.gtf

# Copy local gene annotation (testing only)
# COPY gencode.v40.chr_patch_hapl_scaff.annotation.gtf /home/ubuntu
# COPY gencode.v40.chr_patch_hapl_scaff.annotation.db /home/ubuntu


RUN chmod 644 /home/ubuntu/*

WORKDIR /tools
USER ubuntu
ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]