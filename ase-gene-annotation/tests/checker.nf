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

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
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

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.input_file = ""
params.vcf_file = ""
params.gtf_file = "/home/ubuntu/gencode.v40.chr_patch_hapl_scaff.annotation.gtf"
params.assembly = "GRCh38"
params.expected_output = ""
params.expected_table = ""

include { aseGeneAnnotation } from '../main'


process file_smart_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    path gene_table
    path hap_table
    path expected_output
    path expected_table

  output:
    stdout()

  script:
    """
    diff ${gene_table} ${expected_output}  && ( echo "Test PASSED" ) || ( echo "Test FAILED, output file mismatch." && exit 1 )
    diff ${hap_table} ${expected_table}  && ( echo "Test PASSED" ) || ( echo "Test FAILED, output file mismatch." && exit 1 )
    """
}

workflow checker {
  take:
    input_file
    vcf_file
    expected_output
    expected_table

  main:
    aseGeneAnnotation(
      input_file,
      vcf_file
    )

    file_smart_diff(
      aseGeneAnnotation.out.gene_table,
      aseGeneAnnotation.out.hap_table,
      expected_output,
      expected_table
    )
}


workflow {
  checker(
    file(params.input_file),
    file(params.vcf_file),
    file(params.expected_output),    
    file(params.expected_table)
  )
}
