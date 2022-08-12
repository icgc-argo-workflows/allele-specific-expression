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

import vcf
import argparse
import pandas as pd
import numpy as np
import scipy as sp
import logging
import statsmodels.api as sm

def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%a, %d %b %Y %H:%M:%S', 
        handlers=[
            logging.FileHandler("hap_table.log"),
            logging.StreamHandler()
        ], 
        force=True)

    parser = argparse.ArgumentParser(description='Construct a table counting expression of either parental haplotype across gene-mapped positions.')
    parser.add_argument("-I", "--input_file", required=True, help="Gene-annotated ASE read counter result file.")
    parser.add_argument("-V", "--variants_file", required=True, help="Variants Calling File with genotype annotation (GT).")
    parser.add_argument("-O", "--output_file", required=True, help="Gene table output file.")
    args = parser.parse_args()

    ase = pd.read_csv(args.input_file, sep="\t")
    vcf_reader = vcf.Reader(open(args.variants_file))

    ase_filered = filter_unmatched(ase)
    if len(ase_filered) <= 0:
        logging.warning("No ASE positions. Skipping hap_table process.")
        return

    gts = get_genotypes(vcf_reader)
    genotyped = pd.merge(ase_filered, gts, how='inner', left_on=["position", "contig"], right_on=["position", "contig"])
    gen_filtered = filter_genotypes(genotyped)
    if len(gen_filtered) <= 0:
        logging.warning("No genotype positions. Skipping hap_table process.")
        return

    simplified = gen_filtered.copy().drop(columns=["position", "variantID", "refAllele", "altAllele", "ase_ratio", "ref_bias", "AEI_pval", "AEI_padj"])
    counted = count_hap(simplified)
    grouped = counted.groupby(["gene_id"], as_index=False)
    gene_table = pd.merge(grouped.sum(), grouped.size()).rename(columns={"size": "positions"})
    haplotype_imbalance(gene_table, counted)

    logging.info(f"Exporting data to file {args.output_file}.")
    gene_table.round(4).to_csv(args.output_file, sep="\t", index=False)
    logging.info("Done.")


def filter_unmatched(ase):
    ase_filter = ase["gene_id"].notna()
    logging.info(f"Removing {len(ase) - ase_filter.sum()}/{len(ase)} unmatched positions.")
    return ase[ase_filter]


def filter_genotypes(genotyped):
    phased_filter = genotyped["GT"].str.contains("|")
    print(f"Removing {len(genotyped) - phased_filter.sum()}/{len(genotyped)} unphased positions.")
    genotyped = genotyped[phased_filter]
    dip_het_filter = (genotyped["GT"] == "0|1") | (genotyped["GT"] == "1|0") 
    print(f"Removing {len(genotyped) - dip_het_filter.sum()}/{len(genotyped)} non-diploid positions.")
    return genotyped


def get_genotypes(vcf_reader):
    logging.info(f"Reeading GT info.")
    gts_list = []
    for record in vcf_reader:
        hets = record.get_hets()
        if len(hets) > 0 and hets[0].phased: 
            gts_list.append([record.POS, record.CHROM, hets[0].data.GT])      
    return pd.DataFrame(gts_list, columns=["position", "contig", "GT"])


def count_hap(hap_df):
    hap_df.loc[hap_df["GT"] == "0|1", "hap1"] = hap_df["refCount"] 
    hap_df.loc[hap_df["GT"] == "0|1", "hap2"] = hap_df["altCount"] 
    hap_df.loc[hap_df["GT"] == "1|0", "hap1"] = hap_df["altCount"] 
    hap_df.loc[hap_df["GT"] == "1|0", "hap2"] = hap_df["refCount"] 
    hap_df = hap_df.astype({"hap1": int, "hap2": int})
    return hap_df


def haplotype_imbalance(gene_table, haplotypes):  
    logging.info("Computing the imbalance between haplotypes.")
    gene_table["HSE_ratio"] = haplotypes["hap1"] / haplotypes["totalCount"] 
    gene_table["HEI_pval"] = np.array([sp.stats.binomtest(k=gene_table.hap1[idx], n=gene_table.totalCount[idx], p=.5).pvalue for idx in gene_table.index], dtype="float")
    gene_table["HEI_padj"] = sm.stats.multipletests(gene_table['HEI_pval'], method="fdr_bh")[1]


if __name__ == '__main__':
    main()
