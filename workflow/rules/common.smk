import sys
from snakemake.utils import logger, validate
import itertools
import pandas as pd
import numpy as np

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_PREFIX = Path("https://github.com/snakemake/snakemake-wrappers/raw")
WRAPPER_PREFIX = Path(workflow.wrapper_prefix)
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = Path("https://raw.githubusercontent.com/percyfal/assemblyeval-smk/main/workflow")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

assemblies = pd.read_csv(config["assemblies"], sep="\t").set_index(["species", "version"], drop=False)
assemblies = assemblies.replace({np.nan: None})
assemblies.index.names = ["species", "version"]
validate(assemblies, schema="../schemas/assemblies.schema.yaml")

transcripts = None
if config["transcripts"]:
    transcripts = pd.read_csv(config["transcripts"], sep="\t").set_index("id", drop=False)
    transcripts = transcripts.replace({np.nan: None})
    transcripts.index.names = ["id"]
    validate(transcripts, schema="../schemas/transcripts.schema.yaml")

reads = None
if config["reads"]:
    reads = pd.read_csv(config["reads"], sep="\t").set_index("id", drop=False)
    reads = reads.replace({np.nan: None})
    reads.index.names = ["id"]
    validate(reads, schema="../schemas/reads.schema.yaml")

## Validate genecovr files
for v in config["genecovr"].keys():
    if not v.startswith("dataset"):
        continue
    csvfile = config["genecovr"][v]['csvfile']
    if csvfile is None:
        continue
    data = pd.read_csv(csvfile, header=None).set_index(0, drop=False)
    data.columns = ["dataset", "psl", "assembly", "trxset"]
    assert data["assembly"].str.replace(r".fai", "").isin(assemblies["fasta"]).all(),\
        "some values in 'assembly' column not present in assemblies input file"
    assert data["trxset"].str.replace(r".fai", "").isin(transcripts["fasta"]).all(),\
        "some values in 'trxset' column not present in transcripts input file"
    validate(data, schema="../schemas/genecovr_csv.schema.yaml")

##############################
## Config checks; function defs
##############################
def check_blobdir_keys():
    blobdir_assemblies = []
    for blobdir in config["btk"].keys():
        if not blobdir.startswith("blobdir_"):
            continue
        species, version = re.sub("blobdir_", "", blobdir).split("_")
        if (species, version) not in assemblies.index:
            logger.error(f"error in blobdir configuration: {species}_{version} not in {config['assemblies']}")
            sys.exit(1)




##############################
## Config checks
##############################
check_blobdir_keys()

##############################
## Paths
##############################
__EXTERNAL__ = Path(config["fs"]["external"])
__INTERIM__ = Path(config["fs"]["interim"])
__METADATA__ = Path(config["fs"]["metadata"])
__RAW__ = Path(config["fs"]["raw"])
__REPORTS__ = Path("reports")
__RESOURCES__ = Path(config["fs"]["resources"])
__RESULTS__ = Path("results")

##############################
## Wildcard constraints
##############################
wildcard_constraints:
    external = str(__EXTERNAL__),
    interim = str(__INTERIM__),
    metadata = str(__METADATA__),
    raw = str(__RAW__),
    reports = str(__REPORTS__),
    resources = str(__RESOURCES__),
    results = str(__RESULTS__)

## Assemblies etc
wildcard_constraints:
    assembly = "|".join(f"{species}_{version}" for species, version in assemblies.index),
    blobdir = "[^/]+",
    region = "\d+"

## File extensions
wildcard_constraints:
    fa = "(.fa|.fasta)",
    gz = "(|.gz)"


##################################################
## Functions and utilities
##################################################
def get_assembly(wildcards):
    # In principle could use ensembl wrapper if db not present:
    # file:// would be better
    species, version = wildcards.assembly.split("_")
    return assemblies.loc[species, version]["fasta"]

def get_transcriptome(wildcards):
    value = transcripts.loc[wildcards.transcriptome]["fasta"]
    if isinstance(value, str):
        value = [value]
    return value


##############################
## Pseudo-rule targets
##############################
def get_btk_all(wildcards):
    retval = []
    for v in config["btk"].keys():
        if not v.startswith("blobdir"):
            continue
        if len(config["btk"][v]["fasta"]) > 0:
            retval.append(str(__INTERIM__ / f"btk/{v}/gc.json"))
        for bam in config["btk"][v]["bam"]:
            retval.append(str(__INTERIM__ / f"btk/{v}/{Path(bam).name}.sort_cov.json"))
        if len(config["btk"][v]["bls"]) > 0:
            retval.append(str(__INTERIM__ / f"btk/{v}/bestsumorder_class_score.json"))
    return retval[0]


def get_genecovr_all(wildcards):
    dataset = [x for x in config["genecovr"].keys() if x.startswith("dataset")]
    return expand(f"{str(__RESULTS__)}/genecovr/{{dataset}}/psldata.csv.gz", dataset=dataset)
