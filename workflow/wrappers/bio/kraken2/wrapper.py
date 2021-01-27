__author__ = "Per Unneberg"
__copyright__ = "Copyright 2021, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

options = snakemake.params.get("options", "")
db = snakemake.params.db

seqinput = snakemake.input[0]


shell(
    "kraken2 --db {db} --threads {snakemake.threads} --output {snakemake.output.output} "
    "--report {snakemake.output.report} --report-zero-counts --unclassified-out {snakemake.output.unclassified} "
    "{options} {seqinput} {log}"
)
