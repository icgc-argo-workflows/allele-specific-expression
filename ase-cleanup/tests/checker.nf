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
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/allele-specific-expression.ase-cleanup'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.input_file = ""
params.mapp_file = "/home/ubuntu/k50.umap.bedgraph.gz"
params.genome_file = "/home/ubuntu/GRCh38_Verily_v1.genome.fa.gz.genome"
params.min_mappability = 0.05
params.min_SNP_depth = 16

include { aseCleanup } from '../main'

process file_smart_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    path output_file
    path expected_file
    path log_file
    path vaf_file

  output:
    stdout()

  script:
    """
    diff ${output_file[0]} ${expected_file} && ( echo "Output file matched" ) || ( echo "Test FAILED, output file mismatch." && exit 1 )
    test -e ${log_file} || ( echo "Test FAILED, log missing" && exit 2 )
    test -e ${vaf_file} || ( echo "Test FAILED, VAF png missing" && exit 3 )
    """
}


workflow checker {
  take:
    input_file
    expected_output

  main:
    aseCleanup(
      input_file        
    )

    file_smart_diff(
      aseCleanup.out.output_file,
      expected_output,
      aseCleanup.out.log_file,
      aseCleanup.out.vaf_file
    )
}


workflow {
  checker(
    file(params.input_file),
    file(params.expected_output)
  )
}
