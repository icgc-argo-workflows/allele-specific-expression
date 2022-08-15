# Package ase-gene-annotation


Annotates the ASE gene expression data with the gene labels. 
1. The positions matching introns/exons are annotated with the associated genes.
2. For phased data, a table with haplotype specific expression is created. 

The gene annotation adds following columns:
* `gene_id`: the ENSAMBL gene id,
* `feature`: one of (intron, exon).

The haplotype specific expression table has the following columns:
* `gene_id`: the ENSAMBL gene id,
* `refCount`: the number of reads supporting the reference haplotype,
* `altCount`: the number of reads supporting the alternative haplotype,
* `totalCount`: the total number of reads supporting the haplotype,
* `hap1`: the count of reads from the haplotype 1,
* `hap2`: the count of reads from the haplotype 2,
* `positions`: the number of positions in the gene,
* `HSE_ratio`: the ratio of the haplotype specific expression,
* `HEI_pval`: the p-value of the haplotype expression imballance score,
* `HEI_padj`: the adjusted p-value of the HEI score.

## Package development

The initial version of this package was created by the WorkFlow Package Manager CLI tool, please refer to
the [documentation](https://wfpm.readthedocs.io) for details on the development procedure including
versioning, updating, CI testing and releasing.


## Inputs

* `sample_name.clean`: The result of the ase-cleanup.
* `variants.vcf`: The VCF file with phased data. 


## Outputs

* `sample_name.tsv`: The annotated table.
* `sample_name.gene.log`: Log of the annotation process.
* `sample_name.hap.tsv`: The haplotype specific expression table.
* `sample_name.hap.log`: Log of the haplotype specific expression process.

## Usage

### Run the package directly

With inputs prepared, you should be able to run the package directly using the following command.
Please replace the params file with a real one (with all required parameters and input files). Example
params file(s) can be found in the `tests` folder.

```
nextflow run icgc-argo-workflows/allele-specific-expression/ase-gene-annotation/main.nf -r ase-gene-annotation.v0.1.0 -params-file <your-params-json-file>
```

### Import the package as a dependency

To import this package into another package as a dependency, please follow these steps at the
importing package side:

1. add this package's URI `github.com/icgc-argo-workflows/allele-specific-expression/ase-gene-annotation@0.1.0` in the `dependencies` list of the `pkg.json` file
2. run `wfpm install` to install the dependency
3. add the `include` statement in the main Nextflow script to import the dependent package from this path: `./wfpr_modules/github.com/icgc-argo-workflows/allele-specific-expression/ase-gene-annotation@0.1.0/main.nf`
