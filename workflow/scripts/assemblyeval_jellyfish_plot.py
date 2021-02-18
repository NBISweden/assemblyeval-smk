#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import shutil
import seaborn as sns
import pandas as pd
from snakemake.shell import shell
from snakemake.utils import logger

asm_file = snakemake.input.assembly
read_file = snakemake.input.reads
tsv_file = snakemake.output.tsv

exe = "kmer_count_pairs"
if shutil.which(exe) is None:
    logger.error(
        "[assemblyeval_jellyfish_plot]: no such binary '{exe}'; install in PATH"
    )
    sys.exit(1)

# FIXME: kmer_count_pairs must be built automatically
shell("{exe} {asm_file} {read_file} > {tsv_file}")

with open(tsv_file, "r") as fh:
    df = pd.read_csv(fh, sep="\t")

df.columns = ["asm_mers", "read_mers", "count"]
cmax1 = min(max(df["asm_mers"]), 6)
cmax2 = max(df["read_mers"])
df["asm_mers"] = df["asm_mers"].apply(lambda x: x if x < cmax1 else cmax1)

llab = [f"{x}X" for x in set(df["asm_mers"])]
llab[-1] = f"{llab[-1]}+"

p = sns.lineplot(data=df, x="read_mers", y="count", hue="asm_mers")
p.set(xlim=(0, None), title="k-mer comparison plot")
p.set_xlabel("read k-mer multiplicity")
p.set_ylabel("Number of distinct k-mers")
p.legend(llab, title="Assembly k-mer multiplicity")
fig = p.get_figure()
fig.savefig(snakemake.output.png)
