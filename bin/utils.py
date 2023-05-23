# Python utility functions

import numpy as np
import pandas as pd


def read_maf(file: str) -> pd.DataFrame:
    # count comment lines
    c = 0
    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                c = c + 1

    maf = pd.read_csv(file, header=c, sep="\t", low_memory=False)
    return maf


def writeheader(file: str, outfile: str):
    """Write header from input maf file into outfile."""
    out = open(outfile, "w")
    with open(file, "r") as f:
        for line in f:
            if line.startswith("#"):
                out.write(line)

    out.close()
