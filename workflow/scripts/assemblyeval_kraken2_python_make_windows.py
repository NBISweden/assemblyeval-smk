#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from pybedtools import BedTool


partition = int(snakemake.wildcards.partition)
bed = BedTool(snakemake.input.bed)
out = snakemake.output[0]
nwin = len(bed)
npart = int(nwin / int(snakemake.params.npart))
[imin, imax] = list(range(0, nwin, npart))[partition:(partition+2)]

a = BedTool(bed[imin:imax])
a = a.sequence(fi=snakemake.input.fasta, fo=out)
