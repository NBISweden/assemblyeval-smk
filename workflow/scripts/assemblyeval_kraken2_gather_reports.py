#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Read and merge kraken2 results output
#
from tqdm import tqdm
from functools import reduce
import pandas as pd

def load(fn):
    df = pd.read_table(
        fn, header=None,
        names=['percent', 'total', 'assigned', 'rank', 'taxid', 'name'])
    df.set_index("taxid", drop=False, inplace=True)
    return df

def merge(data, df):
    data.total = data.total + df.total
    data.assigned = data.assigned + df.assigned
    data.percent = data.total / data.total[1] * 100
    return data


reports = snakemake.input.txt

dflist = []
data = load(reports[0])
for fn in tqdm(reports[1:len(reports)]):
    data = merge(data, load(fn))

data.to_csv(snakemake.output.txt, compression="gzip", index=False, sep="\t", header=False, float_format="%6.2f")
