#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
from multiprocessing import Pool
import pybedtools
from tqdm import tqdm
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

faidx = snakemake.input.faidx

regions = []
with open(faidx) as fh:
    for l in tqdm(fh.readlines()):
        fields = l.split("\t")
        i = pybedtools.Interval(fields[0], 0, int(fields[1]))
        regions.append(i)

npart = len(snakemake.output)
try:
    assert len(regions) >= npart
except AssertionError as e:
    logger.warning(
        (
            f"Number of regions smaller than number of partitions: '{len(regions)} < {npart}': "
            f"lower the number of partitions "
        )
    )
    raise

regions.sort(key=lambda r: len(r), reverse=True)
out = [[regions[i]] for i in range(npart)]
outlen = [len(r[0]) for r in out]
for j in tqdm(range(npart, len(regions))):
    imin = outlen.index(min(outlen))
    out[imin].append(regions[j])
    outlen[imin] += len(regions[j])


def _getfasta(args):
    i, out, output, inseq = args
    logger.info(f"saving chunk {i} to {output}...")
    bed = pybedtools.BedTool(out)
    logger.info(f"reading input {inseq}")
    bed.sequence(fi=inseq).save_seqs(output)

logger.info(f"Saving {len(out)} chunks...")
args = [(i, out[i], snakemake.output[i], snakemake.input.seq) for i in range(len(out))]
with Pool(snakemake.threads) as p:
    list(p.imap(_getfasta, args))
