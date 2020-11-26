#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
import gzip
from snakemake.shell import shell
from snakemake.utils import logger

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

seqin = snakemake.input.seq
out = snakemake.output.jf
kmer = snakemake.wildcards.kmer
options = snakemake.params.options
if options == "":
    options = "-s 100M -C"

cat = "cat"

if not isinstance(seqin, list):
    seqin = [seqin]

nzip = sum([re.search(r".gz$", x) is not None for x in seqin])
assert nzip == 0 or len(seqin) == nzip, "infiles must all be unzipped or gzipped!"

if nzip == len(seqin):
    cat = "zcat"

shell(
    "{cat} {seqin} | jellyfish count {options} -m {kmer} -t {snakemake.threads} /dev/fd/0 -o {out} {log}"
)
