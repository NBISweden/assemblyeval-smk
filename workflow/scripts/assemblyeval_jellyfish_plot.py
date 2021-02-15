#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import dna_jellyfish as jf
import numpy as np
import seaborn as sns
from tqdm import tqdm
import pandas as pd

assembly = snakemake.input.assembly
reads = snakemake.input.reads


def load_kmers(fn):
    data = {}
    cmax = 0
    mf = jf.ReadMerFile(fn)
    for mer, count in tqdm(mf):
        data[str(mer)] = count
        cmax = max(cmax, count)
    return cmax, data


cmax1, read_data = load_kmers(reads)
cmax2, assembly_data = load_kmers(assembly)
cmax2 = min(cmax2, 6)

m = np.zeros((cmax1 + 1, cmax2 + 1))

for mer, count in tqdm(assembly_data.items()):
    if count >= m.shape[1]:
        count = m.shape[1] - 1
    if str(mer) in read_data.keys():
        readcount = read_data[str(mer)]
        m[readcount][count] = m[readcount][count] + 1
    else:
        m[0][count] = m[0][count] + 1

for mer, count in tqdm(read_data.items()):
    if str(mer) not in assembly_data.items():
        m[count][0] = m[count][0] + 1

df = pd.DataFrame(m).stack().to_frame().reset_index()
df.columns = ["kcount", "kdistinct", "val"]
cmax = max(df["kdistinct"])
df["kdistinct"] = df["kdistinct"].astype("str").replace({f"{cmax}": f"{cmax}+"})

p = sns.displot(df, x="kcount", weights="val", hue="kdistinct", kind="kde")
p.set(xlim=(0, None), title="k-mer comparison plot")
p.set_axis_labels("k-mer multiplicity", "Number of distinct k-mers")

p.savefig(snakemake.output.png)
df.to_csv(snakemake.output.tsv, index=False, sep="\t")
