#!/usr/bin/env python3
import itertools

import pandas as pd
import snakemake

input = snakemake.input
output = snakemake.output
wildcards = snakemake.wildcards
config = snakemake.config


def genecovr_make_csv_inputfile_dataframe(d, wildcards):
    assembly_keys = config["genecovr"][wildcards.dataset]["assemblies"]
    trx_keys = config["genecovr"][wildcards.dataset]["transcripts"]
    df = pd.concat(
        [
            pd.Series(
                [f"{a}/{b}" for a, b in itertools.product(assembly_keys, trx_keys)]
            ),
            pd.Series(list(d["psl"])),
            pd.Series(list(d["assembly"])),
            pd.Series(list(d["trxset"])),
        ],
        axis=1,
    )
    return df


df = genecovr_make_csv_inputfile_dataframe(dict(input), wildcards)
df.to_csv(output.csv, index=False, header=False)
