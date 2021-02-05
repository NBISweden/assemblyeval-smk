__author__ = "Per Unneberg"
__copyright__ = "Copyright 2021, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

outdir = os.path.dirname(snakemake.output.align)

shell(
    "RepeatMasker -dir {outdir} -a -pa {snakemake.threads} {snakemake.params.options} {snakemake.input} {log} "
)
