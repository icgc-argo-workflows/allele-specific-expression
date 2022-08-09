# Package ase-cleanup

This package runs an ase-cleanup script.

The script takes a GATK ASEReadCounter output and filters reads that do not match the following conditions:
* mappability of a position is lower than `min_mappability` (default: `0.05`)
* the number of reads for a position is lower than `min_SNP_depth` (default: `16`)

Using the above, it produces a tab separated document detailing the results of the ASE analysis with the following result columns:
    * `ase_ratio`: the RAF adjusted for mean bias towards reference
    * `ref_bias`: the ration of reference counts vs total read counts for the particular base pair
    * `AEI_pval`: the resulting p-value of binomial statistical test
    * `AEI_padj`: the p-value corrected using Benjamini/Hochberg false discovery rate correction. The AE is present if `p < 0.5`. 

## Package development

The initial version of this package was created by the WorkFlow Package Manager CLI tool, please refer to
the [documentation](https://wfpm.readthedocs.io) for details on the development procedure including
versioning, updating, CI testing and releasing.


## Inputs

Please list all input parameters


## Outputs

Please list all outputs


## Usage

### Run the package directly

With inputs prepared, you should be able to run the package directly using the following command.
Please replace the params file with a real one (with all required parameters and input files). Example
params file(s) can be found in the `tests` folder.

```
nextflow run icgc-argo-workflows/allele-specific-expression/ase-cleanup/main.nf -r ase-cleanup.v0.1.0 -params-file <your-params-json-file>
```

### Import the package as a dependency

To import this package into another package as a dependency, please follow these steps at the
importing package side:

1. add this package's URI `github.com/icgc-argo-workflows/allele-specific-expression/ase-cleanup@0.1.0` in the `dependencies` list of the `pkg.json` file
2. run `wfpm install` to install the dependency
3. add the `include` statement in the main Nextflow script to import the dependent package from this path: `./wfpr_modules/github.com/icgc-argo-workflows/allele-specific-expression/ase-cleanup@0.1.0/main.nf`
