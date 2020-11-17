__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

assembly = snakemake.input.seq
output = snakemake.output.txt
outdir = os.path.dirname(output)

shell("quast.py {assembly} -o {outdir} {log}")
