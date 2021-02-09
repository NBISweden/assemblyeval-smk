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

output = re.sub(r".gz$", "", snakemake.output.output)
unclassified = re.sub(r".gz$", "", snakemake.output.unclassified)

shell(
    "kraken2 --db {db} --threads {snakemake.threads} --output {output} "
    "--report {snakemake.output.report} --use-names --report-zero-counts --unclassified-out {unclassified} "
    "{options} {seqinput} {log}"
)
shell("gzip -v {output} {log}")
shell("gzip -v {unclassified} {log}")
