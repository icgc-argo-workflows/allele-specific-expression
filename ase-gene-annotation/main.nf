#!/usr/bin/env nextflow

/*
  Copyright (c) 2021, Max Delbrück Center for Molecular Medicine in the Helmholtz Association

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  Authors:
    Adam Streck
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.1.0'  // package version

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/allele-specific-expression.ase-gene-annotation'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific parmas go here, add / change as needed
params.input_file = ""
params.vcf_file = ""
params.gtf_file = "/home/ubuntu/gencode.v40.chr_patch_hapl_scaff.annotation.gtf"
params.assembly = "GRCh38"


process aseGeneAnnotation {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  // input, make update as needed
    path input_file
    path vcf_file


  output:  // output, make update as needed
    path("${input_file.baseName}.tsv"), emit: gene_table
    path("${input_file.baseName}.gene.log"), emit: gene_log
    path("${input_file.baseName}.hap.tsv"), emit: hap_table
    path("${input_file.baseName}.hap.log"), emit: hap_log

  script:
    // add and initialize variables here as needed

    """

    gene_annotation.py -I $input_file -O ${input_file.baseName}.tsv --gtf $params.gtf_file --ref $params.assembly
    mv gene_annotation.log ${input_file.baseName}.gene.log
    hap_table.py -I ${input_file.baseName}.tsv -V $vcf_file -O ${input_file.baseName}.hap.tsv
    mv hap_table.log ${input_file.baseName}.hap.log
    """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  aseGeneAnnotation(
    file(params.input_file),
    file(params.vcf_file)
  )
}
