#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Read and merge kraken2 results output
#
import gzip
import pandas as pd

def merge(data, df):
    data.total = data.total + df.total
    data.assigned = data.assigned + data.assigned

def load(fn):
    return pd.read_table(
        fn, header=None,
        names=['percent', 'total', 'assigned', 'rank', 'taxid', 'name'])

reports = snakemake.input.txt

data = load(reports[0])
for i in range(1, len(reports)):
    df = load(reports[i])
    data = merge(data, df)

data.percent = data.total / data.total[1] * 100

data.to_pickle(snakemake.output.txt)
