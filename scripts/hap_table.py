#!/usr/bin/env python

import argparse
import gzip
from pathlib import Path
import pandas as pd
import numpy as np
import scipy as sp
import logging
import statsmodels.api as sm


OUTPUT_COLUMNS = [
    "gene_id", "phase_set", "hap1", "hap2", "totalCount", "positions",
    "HSE_ratio", "HEI_pval", "HEI_padj"]

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
    parser.add_argument("--sample", help="VCF sample to use. Required for multi-sample VCFs unless --bam-name matches a sample.")
    parser.add_argument("--sample-file", help="File containing the preflight-selected VCF sample.")
    parser.add_argument("--bam-name", help="BAM-derived sample name used to auto-select a matching VCF sample.")
    args = parser.parse_args()

    requested_sample = args.sample
    if args.sample_file:
        selected_sample = Path(args.sample_file).read_text(encoding="utf-8").strip()
        if selected_sample:
            requested_sample = selected_sample

    ase = pd.read_csv(args.input_file, sep="\t")

    ase_filered = filter_unmatched(ase)
    if len(ase_filered) <= 0:
        logging.warning("No gene-matched ASE positions. Writing an empty haplotype table.")
        write_empty_table(args.output_file)
        return

    gts = get_genotypes(args.variants_file, requested_sample, args.bam_name)
    genotyped = pd.merge(ase_filered, gts, how='inner', left_on=["position", "contig"], right_on=["position", "contig"])
    gen_filtered = filter_genotypes(genotyped)
    if len(gen_filtered) <= 0:
        logging.warning("No phased heterozygous genotype positions. Writing an empty haplotype table.")
        write_empty_table(args.output_file)
        return

    simplified = gen_filtered[[
        "gene_id", "position", "refCount", "altCount", "totalCount", "GT", "PS"
    ]].copy()
    counted = count_hap(simplified)
    gene_table = summarize_haplotypes(counted)
    haplotype_imbalance(gene_table)

    logging.info(f"Exporting data to file {args.output_file}.")
    gene_table.round(4).to_csv(args.output_file, sep="\t", index=False)
    logging.info("Done.")


def filter_unmatched(ase):
    ase_filter = ase["gene_id"].notna()
    logging.info(f"Removing {len(ase) - ase_filter.sum()}/{len(ase)} unmatched positions.")
    return ase[ase_filter]


def filter_genotypes(genotyped):
    phased_filter = genotyped["GT"].str.contains("|", regex=False, na=False)
    logging.info("Removing %s/%s unphased positions.",
                 len(genotyped) - phased_filter.sum(), len(genotyped))
    genotyped = genotyped[phased_filter]
    dip_het_filter = genotyped["GT"].isin(["0|1", "1|0"])
    logging.info("Removing %s/%s non-diploid or non-heterozygous positions.",
                 len(genotyped) - dip_het_filter.sum(), len(genotyped))
    return genotyped[dip_het_filter].copy()


def get_genotypes(vcf_path, requested_sample=None, bam_name=None):
    logging.info("Reading GT and phase-set information.")
    gts_list = []
    open_vcf = gzip.open if str(vcf_path).endswith(".gz") else open
    sample_column = None

    with open_vcf(vcf_path, "rt") as vcf_file:
        for line in vcf_file:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                selected_sample = select_sample(samples, requested_sample, bam_name)
                sample_column = header.index(selected_sample)
                logging.info("Using VCF sample %s.", selected_sample)
                continue
            if line.startswith("#"):
                continue
            if sample_column is None:
                raise ValueError("VCF is missing a #CHROM header with sample columns.")

            fields = line.rstrip("\n").split("\t")
            format_keys = fields[8].split(":")
            if "GT" not in format_keys:
                continue
            sample_values = fields[sample_column].split(":")
            gt_index = format_keys.index("GT")
            if gt_index >= len(sample_values):
                continue
            phase_set = "."
            if "PS" in format_keys:
                ps_index = format_keys.index("PS")
                if ps_index < len(sample_values) and sample_values[ps_index] not in {"", "."}:
                    phase_set = sample_values[ps_index]
            gts_list.append([
                int(fields[1]), fields[0], sample_values[gt_index], phase_set
            ])

    return pd.DataFrame(gts_list, columns=["position", "contig", "GT", "PS"])


def select_sample(samples, requested_sample=None, bam_name=None):
    if requested_sample:
        if requested_sample not in samples:
            raise ValueError(
                f"Requested VCF sample '{requested_sample}' is not present. "
                f"Available samples: {', '.join(samples)}")
        return requested_sample
    if bam_name and bam_name in samples:
        return bam_name
    if len(samples) == 1:
        return samples[0]
    if not samples:
        raise ValueError("VCF contains no sample columns.")
    raise ValueError(
        "VCF contains multiple samples and none matches the BAM name. "
        "Specify the sample with --vcf_sample. Available samples: "
        + ", ".join(samples))


def count_hap(hap_df):
    hap_df.loc[hap_df["GT"] == "0|1", "hap1"] = hap_df["refCount"] 
    hap_df.loc[hap_df["GT"] == "0|1", "hap2"] = hap_df["altCount"] 
    hap_df.loc[hap_df["GT"] == "1|0", "hap1"] = hap_df["altCount"] 
    hap_df.loc[hap_df["GT"] == "1|0", "hap2"] = hap_df["refCount"] 
    hap_df = hap_df.astype({"hap1": int, "hap2": int})
    return hap_df


def summarize_haplotypes(counted):
    phase_sets_per_gene = counted.groupby("gene_id")["PS"].nunique()
    split_genes = phase_sets_per_gene[phase_sets_per_gene > 1]
    if not split_genes.empty:
        logging.warning(
            "%s genes span multiple phase sets and will produce one row per phase set.",
            len(split_genes))

    return counted.groupby(["gene_id", "PS"], as_index=False).agg(
        hap1=("hap1", "sum"),
        hap2=("hap2", "sum"),
        totalCount=("totalCount", "sum"),
        positions=("position", "size")).rename(columns={"PS": "phase_set"})


def haplotype_imbalance(gene_table):
    logging.info("Computing the imbalance between haplotypes.")
    gene_table["HSE_ratio"] = gene_table["hap1"] / gene_table["totalCount"]
    gene_table["HEI_pval"] = np.array([sp.stats.binomtest(k=gene_table.hap1[idx], n=gene_table.totalCount[idx], p=.5).pvalue for idx in gene_table.index], dtype="float")
    gene_table["HEI_padj"] = sm.stats.multipletests(gene_table['HEI_pval'], method="fdr_bh")[1]


def write_empty_table(output_file):
    pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(output_file, sep="\t", index=False)


if __name__ == '__main__':
    main()
