#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (c) 2021, ICGC ARGO
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


import pandas as pd
import numpy as np
import argparse
import scipy as sp
import statsmodels.api as sm
import logging
import matplotlib.pyplot as plt


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%a, %d %b %Y %H:%M:%S',
        handlers=[
            logging.FileHandler("ase_cleanup.log"),
            logging.StreamHandler()
        ],
        force=True)

    parser = argparse.ArgumentParser(description="merge mappability table, filter mappability score, calculate allelic imbalance and reference bias, test heterozygosity, and clean up ase table by given threshold")
    parser.add_argument("--ase", dest="ase",
                        help="input ASE table created by ASEReadCounter")
    parser.add_argument("--ref_ratio", dest="ref_ratio", type=float,
                        help="Numeric, ref_ratio for calculating allelic imbalance in binomial test, default=pipeline-calculated ref bias from data")
    parser.add_argument("--min_SNP_depth", dest="filter_total_read", default=20, type=int,
                        help="Numeric, threshold of total read at SNPs")
    parser.add_argument("--mappability", dest="mapp",
                        help="mappability file of the format (contig, pos, mappability)")
    parser.add_argument("--filter_mapp", dest="filter_mapp_score", default=0.05, type=float,
                        help="Numeric, minimum required mappability, only used if mappability file is provided")
    parser.add_argument("--pvalue_het", dest="het_perror", type=float, 
                        help="Numeric, p value threshold for testing heterozygosity, suggest=0.02, default=otherBases.sum()/rawDepth.sum()")
    parser.add_argument("--output", dest="output_file",
						help="output table")
    parser.add_argument("--plot", dest="plot_file",
						help="VAF plot file")
    args = parser.parse_args()

    ase = read_ase(args.ase)
    if args.mapp is not None:
        ase = add_mapp(ase, args.mapp)
    ref_source_cutoff = args.filter_total_read // 2
    ref_bias_table = calc_ref_bias(ase, ref_source_cutoff) if args.ref_ratio is None else insert_ref_bias(args.ref_ratio)
    ase_rb = pd.merge(ase, ref_bias_table, how="left", left_on=["refAllele", "altAllele"], right_on=["refAllele", "altAllele"])
    ase_rb["ref_bias"].fillna(float(0.5), inplace=True)
    ase_rb.set_index(ase.index, inplace=True)
    ase_ai = allelic_imbalance(ase_rb)
    perror = args.het_perror if args.het_perror else ase.otherBases.sum() / ase.rawDepth.sum()
    ase_het = het_test(ase_ai, perror)
    logging.info(f"args.filter_mapp_score = {args.filter_mapp_score}")
    ase_clean = clean_up(ase_het, args.filter_mapp_score, args.filter_total_read, perror)
    ase_clean.to_csv(args.output_file, sep="\t", index=True, header=True)
    vaf_plot(ase_clean, args.plot_file)


def vaf_plot(ase, plot_file):
    step_val = 1/8
    plt.xticks(np.arange(0, 1+step_val, step=step_val))
    plt.hist(ase["ase_ratio"], bins=32)
    plt.xlabel("B-Allele Frequency")
    plt.ylabel("Number of Positions")
    plt.savefig(plot_file)
    

def read_ase(ase_file):
    ase = pd.read_csv(ase_file, sep="\t", index_col=(0, 1))

    if len(ase) <= 0:
        raise Exception("The ASE table is empty.")

    # the reference allele specific expression ratio
    ase = ase[ase["totalCount"] > 0]
    ase['ase_ratio'] = ase['refCount']/ase['totalCount']

    return ase


# Castel et al. 2015: In previous work [5, 8, 29–31] and in this paper, unless mentioned otherwise, we remove about
# 20 % of het-SNPs that either fall within regions of low mappability (ENCODE 50 bp mappability score < 1) or show
# mapping bias in simulations [27]. This reduces the number of sites with strong bias by about 50 %.
def add_mapp(ase, mappability_file):
    mapp = pd.read_csv(mappability_file, sep="\t", index_col=(0, 1))
    ase['mappability'] = pd.to_numeric(mapp['mappability'], errors='coerce')
    return ase


# Castel et al. 2015: he genome-wide reference ratio remaining slightly above 0.5 indicates residual bias (Figure S6a
# in Additional file 6). Using this ratio as a null in statistical tests instead of 0.5 [5, 6] can improve results (
# Figure S6b–e in Additional file 6).
def calc_ref_bias(ase, cutoff):
    logging.info("Computing reference bias")  # equals to mapping bias

    ref_bias_table = ase.query(f"altCount >= {cutoff} and refCount >= {cutoff}")\
        .groupby(['refAllele', 'altAllele'], as_index=False)['ase_ratio'] \
        .mean() \
        .rename(columns={"ase_ratio": "ref_bias"})

    logging.info("Estimated mean reference bias: %.4f" % ref_bias_table["ref_bias"].mean())

    return ref_bias_table


def insert_ref_bias(bias_param):
    bases = pd.DataFrame(list("ACGT"))
    ref_bias_table = pd.merge(bases, bases, how="cross").rename(columns={"0_x": "refAllele", "0_y": "altAllele"})
    ref_bias_table["ref_bias"] = float(bias_param)
    logging.info("Provided mean reference bias: %.4f" % ref_bias_table["ref_bias"].mean())
    return ref_bias_table


def allelic_imbalance(ase):
    logging.info("Computing allelic expression imballance")

    # allelic epxression (effect size)
    ase['AEI_pval'] = np.array([sp.stats.binom_test(ase.refCount[idx], n=ase.totalCount[idx], p=ase.ref_bias[idx]) for idx in ase.index], dtype="float")
    ase['AEI_padj'] = sm.stats.multipletests(ase['AEI_pval'], method="fdr_bh")[1]

    return ase


# Castel et al. 2015: Highly covered sites are rarely strictly monoallelic even in a homozygous state due to rare
# errors in sequencing and alignment (Figure S7b in Additional file 8). Thus, we propose a genotype error filter
# where the average amount of such sequencing noise per sample is first estimated from alleles other than reference (
# REF) or alternative (ALT) (Figure S7c in Additional file 8). Then, binomial testing is used to estimate if the
# counts of REF/ALT alleles are significantly higher than this noise, and sites where homozygosity cannot be thus
# rejected are flagged as possible errors (Fig. 4b). Additionally, it may be desirable to flag fully monoallelic
# sites with low total counts, where homozygosity cannot be significantly rejected, but heterozygosity is not
# supported either. This test can also be applied to study designs with RNA-seq data from multiple samples (e.g.,
# tissues or treatments) of a given individual, genotyped only once, since genotyping error causes consistent
# monoallelic expression in every tissue. In the Geuvadis data set with 1000 Genomes phase 1 genotypes and sites
# covered by eight or more reads, an average of 4.3 % of sites per sample are excluded by these criteria [1 % false
# discovery rate (FDR)].
def het_test(ase, perror):
    het_pvals = [sp.stats.binom_test(np.minimum(ase.altCount[idx], ase.refCount[idx]), n=ase.totalCount[idx], p=perror, alternative="greater") for idx in ase.index]
    het_padj = sm.stats.multipletests(np.array(het_pvals, dtype="float"), method="fdr_bh")[1]

    ase['het_padj'] = het_padj

    return ase


def clean_up(ase, filter_mapp_score, filter_total, perror):
    ai_padj = 0.05
    het_test_fdr = 0.05
    has_mappability = 'mappability' in ase

    is_mapp = ase['mappability'] >= filter_mapp_score if has_mappability else True
    is_enough_reads = ase["totalCount"] >= filter_total
    is_allelic_imbalance = ase["AEI_padj"] < ai_padj
    is_statistically_biallelic = ase['het_padj'] < het_test_fdr

    nsites_original = ase.shape[0]

    if has_mappability:
        logging.info("%d / %d (%.2f%%) of sites removed due to mappability (mappability < %.2f)."
                    % (np.sum(~is_mapp), nsites_original, 100 * np.sum(~is_mapp) / nsites_original, filter_mapp_score))
    logging.info("%d / %d (%.2f%%) of sites removed due to not have enough totalCount (< %d reads)"
                 % (np.sum(~is_enough_reads), nsites_original, 100*np.sum(~is_enough_reads) / nsites_original, filter_total))
    logging.info("%d / %d (%.2f%%) of sites removed due to being staistically homozygous (het_test_fdr > %.2f, perror estimated at %.4f)."
                 % (np.sum(~is_statistically_biallelic), nsites_original, 100 * np.sum(~is_statistically_biallelic) / nsites_original, het_test_fdr, perror))
    logging.info("%d / %d (%.2f%%) of do not have allelic  imbalance (ai_padj > %.2f)."
                 % (np.sum(~is_allelic_imbalance), nsites_original, 100 * np.sum(~is_allelic_imbalance) / nsites_original, ai_padj))

    ase = ase[is_mapp & is_statistically_biallelic & is_enough_reads]
    nsites_now = ase.shape[0]

    logging.info("%d / %d (%.2f%%) of sites removed in total." % (nsites_original -
                 nsites_now, nsites_original, (1-nsites_now/nsites_original) * 100))

    if has_mappability:
        ase = ase.drop(columns=["mappability"])
    return ase.drop(columns=["lowMAPQDepth", "lowBaseQDepth", "rawDepth", "otherBases", "improperPairs", "het_padj"])


if __name__ == "__main__":
    main()
