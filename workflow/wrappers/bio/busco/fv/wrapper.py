__author__ = "Per Unneberg"
__copyright__ = "Copyright 2023, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

from snakemake.shell import shell

shell("busco --download")
