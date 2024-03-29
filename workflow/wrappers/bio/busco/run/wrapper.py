__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

conda_prefix = os.getenv("CONDA_PREFIX")
augustus = os.path.join(conda_prefix, "etc/conda/activate.d/augustus.sh")
if os.path.exists(augustus):
    shell("source {augustus}")

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

options = snakemake.params.get("options", "")
lineage = snakemake.wildcards.lineage
mode = snakemake.wildcards.mode
assembly = snakemake.input.assembly


out_path = os.path.dirname(os.path.dirname(os.path.dirname(snakemake.output.tsv)))
outdir = snakemake.wildcards.mode

buscolineage = re.sub("_odb[0-9]+$", "", lineage)
shell(
    "busco --in {assembly} --lineage {buscolineage} --force "
    "--out {outdir} --out_path {out_path} --cpu {snakemake.threads} --mode {mode} "
    "{options} {log}"
)

short_source = snakemake.output.summary
short_multiqc = os.path.join(
    out_path,
    outdir,
    f"run_{snakemake.wildcards.lineage}",
    f"short_summary_{snakemake.wildcards.assembly}.txt",
)
shell("cp {short_source} {short_multiqc} {log}")
