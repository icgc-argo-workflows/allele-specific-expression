# Pipeline Allele Specific Expression

This is a [Nextflow](https://nextflow.io/) pipeline for analysis of allele specific expression (ASE) based on RNA-seq data. 

## Pipeline

The pipeline processes a reads file (`BAM`/`SAM`) and its accompanying variant call files (`VCF`) together with their matching reference (`FA`) file and produces, for each SNP, the allele-specific expression ratio and the probability that true ASE is occuring.

If a matching genotype (`GTF`) file is provided, the positions are also mapped to genes. 

If the data are additionally phased, the haplotype-specific expression for each gene is computed.

The whole pipeline operates as illustrated:

![Pipeline illustration](./readme/ASE_Pipeline.drawio.png)

The pipeline can be modified using the following quality control parameters:

* QC1 paramters (applied before ASE read counter)
    * base read depth 
    * read mapping quality 
    * base calling quality 
* QC2 parameters (applied after ASE read counter)
    * read mappability (optional)

## Output

For a file `read_name.bam` we obtain the following outputs: 

* `read_name.tsv`: tab separated document detailing the results of the ASE analysis with the following result columns:
    * `ase_ratio`: the RAF adjusted for mean bias towards reference
    * `ref_bias`: the ration of reference counts vs total read counts for the particular base pair
    * `AEI_pval`: the resulting p-value of binomial statistical test
    * `AEI_padj`: the p-value corrected using Benjamini/Hochberg false discovery rate correction. AE is considered significant if `p < 0.05`.
    * [`gene_id`]: if a genome file is provided, maps to an Ensembl gene id
    * [`feature`]: if a genome file is provided, maps to an Ensembl gene feature (exon/intron)
* `read_name.ase.log`: a log file of the ASE calculation and filtering
* `read_name.vaf.png`: a histogram of `ase_ratio` occurrences. In a healthy sample the values should be around 0.5.
* [`read_name.hap.tsv`]: if a genome file is provided and the data are phased, the results of ASE mapped to genes, with the following result columns:
    * `positions`: how many positions are covered by a gene
    * `HSE_ratio`: ratio of the first haplotype vs. total
    * `HEI_pval`: the resulting p-value of binomial statistical test for the gene
    * `HEI_padj`: the FDR B/H p-value correction
* [`read_name.gene_annotation.log`]: a log file of the GTF gene-matching calculation
* [`read_name.hap.log`]: a log file of the haplotype-specific expression calculation

## Requirements

* Java Runtime
* Nextflow version 21.0 or newer
* Conda or Docker

### Nextflow installation

Nextflow can be installed into the repository directly using the command (from the `pipeline-ase` folder):

`wget -qO- https://get.nextflow.io | bash`

### Docker installation

From the root of the repositiory

```
cd pipeline-ase
docker build . -t pipeline-ase
```

### Conda installation

Conda installs itself on the first run of the pipeline. 

## Execution

To execture run 
`nextflow run pipeline-ase.nf [parameters]`

### Output

By default the results are written to the `./out` folder.

### Test run

Random testing data are provided as a part of the package. The default test execution can be done as following:

* Download the reference data using `sh ./DownloadRefData.sh`.
* Execute the pipeline using either Conda or Docker as follows:
    * [Conda] `nextflow -C test_data.config run pipeline-ase.nf -with-conda ase-environment.yaml`
    * [Docker] `nextflow -C test_data.config run pipeline-ase.nf -with-docker pipeline-ase`

**NOTE:** The parameters for the execution are stored in the Nextflow configuration file `test_data.config`.

### Parameters

#### Mandatory
* `--fa path` reference file either as `.fa`  or `.fa.gz`.
* `--bam path` alignment files either as `.bam` or `.sam`. This can be a glob pattern (e.g. `sample*.bam`). All files matching the pattern are used then.

**NOTE:** Read and variants are aligned by name. For example `sample5.bam` is matched to `sample5.vcf`. If you have single vcf for all samples, use the `vcf` parameter detailed below.

#### Default 
* `--min_depth int` minimum read depth for a position, if not satisfied, the position is discarded, `default=16`,
* `--min_mapping_quality int` minimum mapping quality, if not satisfied, the read is discarded, `default=20`,
* `--min_base_quality int` minimum base call quality, if not satisfied, the base is not counted, `default=10`,

#### Optional
* `--vcf path` variant call files either as `.vcf` or `.vcf.gz`. A single file is shared by all BAMs; multiple files from a glob (e.g. `sample*.vcf`) are matched to BAMs by basename. If the path is not provided, the read files are matched as described above.
* `--out_dir path` a path to a folder where the output is stored
* `--assembly name` genome assembly name used by PyEnsembl. Required when `--gtf` is provided (for example, `GRCh37` or `GRCh38`).
* `--gtf path` GTF annotation used for gene matching. `chr1`/`1` style contig aliases are resolved automatically when possible.
* `--vcf_sample name` sample column to use for haplotype calculations with a multi-sample VCF. If omitted, the BAM-derived name is used when it matches a VCF sample; single-sample VCFs are selected automatically.
* `--mapp path` a mappability file in the [bedgraph format](https://genome.ucsc.edu/goldenPath/help/bedgraph.html), which is then used to filter out positions in regions of low quality. Contigs must match the reference
* `--min_mappability float` a number between `[0-1]` describing the minimum mapping quality if the mappability file is present, bases that have lower than specified mappability are removed, `default=0.05`

## Reference downloads

The following human genome files have been tested with the pipeline:

#### reference
`wget https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa`

#### mappability file
`wget https://bismap.hoffmanlab.org/raw/hg38/k50.umap.bedgraph.gz`

#### gene mapping
`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz`

## Contact
Email questions, feature requests and bug reports to **Adam Streck, adam.streck@mdc-berlin.de**.

## License
Pipeline-ASE is available under the MIT License. 

## Please cite
PCAWG Transcriptome Core Group., Calabrese, C., Davidson, N.R. et al. 
Genomic basis for RNA alterations in cancer. Nature 578, 129–136 (2020). 
[https://doi.org/10.1038/s41586-020-1970-0](https://www.nature.com/articles/s41586-020-1970-0)

### Further reading

[Tools and best practices for data processing in allelic expression analysis, Castel et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)
