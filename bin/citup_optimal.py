#!/usr/bin/env python

import sys
import argparse
import pandas as pd


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="h5 file from citup"
    )
    return parser.parse_args()


def main():
    args = _parse_args()
    optimal = pd.read_hdf(args.input, "/results/optimal")
    optimal = optimal.index[optimal == True][0]

    print(optimal)


if __name__ == "__main__":
    main()
