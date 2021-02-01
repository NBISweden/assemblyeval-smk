__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

from os import path

from snakemake.shell import shell


inputs = snakemake.input
output_dir = path.dirname(snakemake.output[0])
output_name = path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
modules = [f"-m {k}" for k, v in inputs.items() if len(v) > 0]

shell(
    "multiqc"
    " {snakemake.params}"
    " --force"
    " -o {output_dir}"
    " -n {output_name}"
    " {modules}"
    " {inputs}"
    " {log}"
)
