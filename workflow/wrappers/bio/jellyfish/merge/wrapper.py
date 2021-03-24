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

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

jfin = snakemake.input.jf
jfout = snakemake.output.jf
kmer = snakemake.wildcards.kmer
options = snakemake.params.options
out = jfout

if not isinstance(jfin, list):
    seqin = [jfin]

npartitions = snakemake.params.npartitions
if npartitions > 1:
    cmd = "jellyfish merge {jfin} -o {out} {log}"
else:
    cmd = "cp {jfin} {jfout}"

shell(cmd)
