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
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/allele-specific-expression.ase-cleanup'
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
params.mapp_file = "/home/ubuntu/k50.umap.bedgraph.gz"
params.genome_file = "/home/ubuntu/GRCh38_Verily_v1.genome.fa.gz.genome"
params.min_mappability = 0.05
params.min_SNP_depth = 16


process aseCleanup {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

    input:
    path(ase)

    output:
    path("*.clean"), emit: output_file
    path("*.log"), emit: log_file
    path("*.vaf.png"), emit: vaf_file

    script:
      var name = ase.baseName
      """ 
      awk -v OFS=\"\\t\" \"{print \\\$1,\\\$2-1,\\\$2 }\" $ase | tail -n+2 | sort -k1,1V -s > ${name}.bed
      echo -e \"contig\tpos\tmappability\" > ${name}.mapp
      bedtools map -a ${name}.bed -b $params.mapp_file -o min -c 4 -g $params.genome_file | cut -f1,3,4 >> ${name}.mapp
      main.py --ase $ase --min_SNP_depth $params.min_SNP_depth  --output ${ase.baseName}.clean --mappability ${name}.mapp --filter_mapp $params.min_mappability --plot ${ase.baseName}.vaf.png
      mv ase_cleanup.log ${ase.baseName}.ase.log
      """
}


// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  aseCleanup(
    file(params.input_file)
  )
}
