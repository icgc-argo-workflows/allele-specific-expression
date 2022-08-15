# Package ase-read-counter


Runs GATK4 ASEReadCounter tool on the provided BAM alignment and the provided germline VCF file. Uses the ICGC AGRO reference RNA genome sequence.


## Package development

The initial version of this package was created by the WorkFlow Package Manager CLI tool, please refer to
the [documentation](https://wfpm.readthedocs.io) for details on the development procedure including
versioning, updating, CI testing and releasing.


## Inputs

* `sample_name.bam`: BAM alignment file
* `mutations.vcf`: SNPs for the given genome.


## Outputs

* `sample_name.read`: The output of the GATK4 ASEReadCounter tool.


## Usage

### Run the package directly

With inputs prepared, you should be able to run the package directly using the following command.
Please replace the params file with a real one (with all required parameters and input files). Example
params file(s) can be found in the `tests` folder.

```
nextflow run icgc-argo-workflows/allele-specific-expression/ase-read-counter/main.nf -r ase-read-counter.v0.1.0 -params-file <your-params-json-file>
```

### Import the package as a dependency

To import this package into another package as a dependency, please follow these steps at the
importing package side:

1. add this package's URI `github.com/icgc-argo-workflows/allele-specific-expression/ase-read-counter@0.1.0` in the `dependencies` list of the `pkg.json` file
2. run `wfpm install` to install the dependency
3. add the `include` statement in the main Nextflow script to import the dependent package from this path: `./wfpr_modules/github.com/icgc-argo-workflows/allele-specific-expression/ase-read-counter@0.1.0/main.nf`
