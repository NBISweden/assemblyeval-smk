#!/usr/bin/env python3
__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

seqin = snakemake.input.seq
jfout = snakemake.output.jf
kmer = snakemake.wildcards.kmer
options = snakemake.params.options
tmpdir = snakemake.params.tmpdir
if tmpdir is not None:
    out = os.path.expandvars(os.path.join(tmpdir, os.path.basename(jfout)))
else:
    out = jfout

if options == "":
    options = "-s 100M -C"

cat = "cat"

if not isinstance(seqin, list):
    seqin = [seqin]

nzip = sum([re.search(r".gz$", x) is not None for x in seqin])
assert nzip == 0 or len(seqin) == nzip, "infiles must all be unzipped or gzipped!"

if nzip == len(seqin):
    cat = "zcat"

cmd = "{cat} {seqin} | jellyfish count {options} -m {kmer} -t {snakemake.threads} /dev/fd/0 -o {out} {log}"

if tmpdir is not None:
    cmd = f"{cmd}; cp {{out}} {{jfout}} {{log}}"

shell(cmd)
