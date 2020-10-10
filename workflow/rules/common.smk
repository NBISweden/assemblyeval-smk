from snakemake.utils import validate
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
    assembly = "|".join(f"{species}_{version}" for species, version in assemblies.index)

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
