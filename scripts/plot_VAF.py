#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description='Crate a Variable Allele Frequency plot over ASEReadCounter output data.')
    parser.add_argument("-I", "--input_file", required=True, help="ASE read counter result file.")
    parser.add_argument("-O", "--output_file", required=True, help="VAF plot file")
    args = parser.parse_args()

    ase = pd.read_csv(args.input_file, sep="\t")
    ase = ase[ase["totalCount"] > 0]

    step_val = 1/8
    plt.xticks(np.arange(0, 1+step_val, step=step_val))
    plt.hist(ase["ase_ratio"], bins=32)
    plt.xlabel("B-Allele Frequency")
    plt.ylabel("Number of Positions")
    plt.savefig(args.output_file)


if __name__ == '__main__':
    main()