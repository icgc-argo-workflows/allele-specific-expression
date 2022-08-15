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

nextflow.enable.dsl = 2
version = '0.1.0'  // package version

// universal params go here, change default value as needed
params.container = ""
params.container_registry = ""
params.container_version = ""
params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir

// tool specific parmas go here, add / change as needed
params.bam = ""
params.vcf = ""
params.cleanup = true

include { aseReadCounter } from './wfpr_modules/github.com/icgc-argo-workflows/allele-specific-expression/ase-read-counter@0.1.0/ase-read-counter'
include { aseCleanup } from './wfpr_modules/github.com/icgc-argo-workflows/allele-specific-expression/ase-cleanup@0.1.0/ase-cleanup'
include { aseGeneAnnotation } from './wfpr_modules/github.com/icgc-argo-workflows/allele-specific-expression/ase-gene-annotation@0.1.0/ase-gene-annotation'


// please update workflow code as needed
workflow Ase {
  take:  // update as needed
    bam
    vcf


  main:  // update as needed
    read_out = aseReadCounter(bam, vcf)
    clean_out = aseCleanup(read_out.output_file)
    annotate_out = aseGeneAnnotation(clean_out.output_file)

  emit:  // update as needed
    output_ase = annotate_out.out.gene_table
    output_hse = annotate_out.out.hap_table
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  Ase(
    file(params.bam),
    file(params.vcf)
  )
}