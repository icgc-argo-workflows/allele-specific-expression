#!/usr/bin/env python

"""
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
"""

import os
import argparse
import logging
import pandas as pd
import pyensembl


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%a, %d %b %Y %H:%M:%S', 
        handlers=[
            logging.FileHandler("gene_annotation.log"),
            logging.StreamHandler()
        ], 
        force=True)

    # Register arguments
    parser = argparse.ArgumentParser(description='Command line utility to '
                                     'create a gene table with contig, pos, gene_id, ... columns from a '
                                     'genomic postion file with contig, pos, ... columns. Genomic positions '
                                     'overlapping multiple genes will be duplicated for each gene. Genomic '
                                     'positions without overlap will not be included in output.')
    parser.add_argument("--gtf", required=True, help="Annotation file as GFT format, in case pyensembl install failed")
    parser.add_argument("-I", "--genomic_position_file", required=True, help="File with single base genomic position coordinates as first as second column. Format (tab-separated, with header): contig, pos, ...")
    parser.add_argument("-O", "--output_file", required=True, help="Gene table output file. Format (tab-separated, with header): contig, pos, gene_id, ...")
    parser.add_argument("--ref", help="Reference name")

    args = parser.parse_args()

    ase = pd.read_csv(args.genomic_position_file, sep="\t")

    indices = ["contig", "position"]

    gpos = ase[indices].copy()
    gpos["contig"] = gpos["contig"].transform(lambda c: c[3:] if c[:3] == "chr" else c)

    data = pyensembl.Genome(reference_name=args.ref, annotation_name='genome_annotation', gtf_path_or_url=args.gtf)
   
    gene_table = create_gene_table(data, ase, gpos)

    result = pd.merge(ase, gene_table, how='left', left_on=indices, right_on=indices)

    ## Export resulting gene table
    logging.info("Exporting data to file %s.", args.output_file)
    result.to_csv(args.output_file, sep="\t", index=False)
    logging.info("Done.")


def create_gene_table(data, ase, gpos):
    biotypes_include = {}
    biotypes_exclude = {}
    columns = ["contig", "position", "gene_id", "feature"]
    gene_table = pd.DataFrame(columns=columns)
    n_sites = len(gpos)

    idx = 0
    stats_sites = 0
    stats_multigene_sites = 0
    stats_duplicates = 0
    stats_sites_lost = 0

    logging.info('Processing sites.')

    for index, row in ase.iterrows():
        contig = row["contig"]
        # contig_str = contig[3:] if contig[:3] == "chr" else contig
        pos = row["position"]
        genes = data.genes_at_locus(contig=contig, position=pos)
        exons = data.exons_at_locus(contig=contig, position=pos)
        exons_gene_ids = set([exon.gene_id for exon in exons])

        start_idx = idx
        for gene in genes:
            feature = 'exon' if gene.gene_id in exons_gene_ids else 'intron'

            # Add gene info to table
            gene_cols = [contig, pos, gene.gene_id, feature]
            gene_table.loc[idx] = gene_cols
            idx += 1

        # Did we add more than one gene at this site?
        if idx - start_idx > 1:
            stats_duplicates += idx - start_idx - 1
            stats_multigene_sites += 1
        # Did we not add any gene at this site?
        elif idx == start_idx:
            stats_sites_lost += 1
        stats_sites += 1
        if(stats_sites % 10000 == 0):
            logging.info('%s sites (%.2f%%) processed.', stats_sites,
                         float(stats_sites)/float(n_sites)*100)

    int_cols = ['position', 'start', 'end']
    cols = list(gene_table.columns.values)
    for col in int_cols:
        if col in cols:
            # Convert POS columns to integer
            gene_table[col] = gene_table[col].astype(int)

    # Log stats
    logging.info('%s sites processed.', stats_sites)
    logging.info('%s multigene sites found.', stats_multigene_sites)
    logging.info('%s position duplicates added at multigene sites.', stats_duplicates)
    logging.info('%s sites lost.', stats_sites_lost)

    return(gene_table)


if __name__ == '__main__':
    main()
