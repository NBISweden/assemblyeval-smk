__author__ = "Per Unneberg"
__copyright__ = "Copyright 2020, Per Unneberg"
__email__ = "per.unneberg@scilifelab.se"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell

# Check genome size
gmap = "gmap"
with open(snakemake.input.log) as fh:
    m = re.search(r"Total genomic length = (\d+) bp", "\n".join(fh.readlines()))
    try:
        if int(m.group(1)) > 2**32:
            gmap = "gmapl"
    except Exception as e:
        print(e)
        raise

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

dbok = snakemake.input.db
transcriptome = snakemake.input.transcriptome
output = snakemake.output.res


def compose_input_gz(filename):
    if filename.endswith(".gz"):
        return f"<(gzip --decompress --stdout {filename})"
    return filename


input_files = [compose_input_gz(fn) for fn in transcriptome]
dirname = os.path.dirname(dbok)
db = f"{snakemake.wildcards.assembly}.db"

shell(
    "{gmap} -t {snakemake.threads} "
    "--dir {dirname} "
    "--db {db} "
    "-f 1 "
    "{input_files} "
    "> {output} "
    "{log}"
)
