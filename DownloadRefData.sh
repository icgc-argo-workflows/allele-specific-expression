#!/usr/bin/env bash

mkdir -p refs
cd refs
wget https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa
bgzip GRCh38_Verily_v1.genome.fa
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
wget https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bedgraph.gz 