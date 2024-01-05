#!/usr/bin/env python3
import pybedtools
from snakemake.utils import logger
from tqdm import tqdm

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

faidx = snakemake.input.faidx
partition = int(snakemake.wildcards.partition)
analysis = snakemake.wildcards.analysis

regions = []
with open(faidx) as fh:
    for line in tqdm(fh.readlines()):
        fields = line.split("\t")
        i = pybedtools.Interval(fields[0], 0, int(fields[1]))
        regions.append(i)

bed = pybedtools.BedTool(regions)
npart = snakemake.params.npartitions

try:
    assert len(bed) >= npart
except AssertionError:
    logger.warning(
        f"Number of regions smaller than number of partitions: '{len(bed)} < {npart}': "
        f"lower the number of partitions "
    )
    raise

try:
    assert partition < npart
except AssertionError:
    logger.error(
        f"partition number {partition} larger than the maximum number of partitions {npart}"
    )
    raise

logger.info(f"saving chunk {partition} to {snakemake.output.fasta}...")
out = bed.splitbed(n=npart).bedtools
bed = out[partition]
logger.info(f"reading input {snakemake.input.seq}")
bed.sequence(fi=snakemake.input.seq).save_seqs(snakemake.output.fasta)
