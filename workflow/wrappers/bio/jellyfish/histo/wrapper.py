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

counts = snakemake.input.counts
hist = snakemake.output.hist
options = snakemake.params.options
tmpdir = snakemake.params.tmpdir

if tmpdir is not None:
    incounts = os.path.join(tmpdir, os.path.basename(counts))
    cmd = "cp {counts} {incounts} {log};"
else:
    incounts = counts
    cmd = ""

cmd = f"{cmd} jellyfish histo {{options}} {{incounts}} -o {{hist}} {{log}}"

shell(cmd)
