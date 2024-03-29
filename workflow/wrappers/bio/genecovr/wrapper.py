__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
from snakemake.shell import shell

options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

analysis = snakemake.wildcards.analysis
outdir = os.path.join("results", "genecovr", analysis)

shell(
    "genecovr -p {snakemake.threads} "
    "{snakemake.params.options} "
    "-d {outdir} "
    "{snakemake.input.csv} "
    "{log}"
)
