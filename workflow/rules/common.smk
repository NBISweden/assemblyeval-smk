import sys
import re
import os
import itertools
import urllib
import pandas as pd
import numpy as np
from snakemake.utils import validate
from snakemake.logging import logger
from snakemake.utils import logger, validate

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_VERSION = "0.67.0"
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
SMK_WRAPPER_PREFIX = f"{SMK_WRAPPER_PREFIX_RAW}/{SMK_WRAPPER_VERSION}"
WRAPPER_PREFIX = workflow.wrapper_prefix
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/percyfal/assemblyeval-smk/main/workflow/wrappers"

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


def _read_tsv(infile, index, schema):
    if infile is None:
        return None
    df = pd.read_csv(infile, sep="\t").set_index(index, drop=False)
    df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schema)
    return df

assemblies = _read_tsv(config["assemblies"], ["species", "version"],
                       "../schemas/assemblies.schema.yaml")
datasources = _read_tsv(config["datasources"], ["data"],
                        "../schemas/datasources.schema.yaml")
transcripts = _read_tsv(config["transcripts"], ["id"],
                        "../schemas/transcripts.schema.yaml")
reads = _read_tsv(config["reads"], ["id"],
                  "../schemas/reads.schema.yaml")


# Save current base dir for later validation in functions
BASEDIR=workflow.current_basedir

##################################################
## Core configuration
##################################################
include: "core/config.smk"

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
    assembly = "|".join(make_assembly_ids()),
    blobdir = "[^/]+",
    region = "[0-9]+"

## File extensions
wildcard_constraints:
    fa = "(.fa|.fasta)",
    gz = "(|.gz)"

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_workdir__"] = os.getcwd()
config["__worfklow_commit__"] = None

try:
    with cd(workflow.basedir, logger):
        logger.info(f"cd to {workflow.basedir}")
        out = sp.check_output(["git", "rev-parse", "--short", "HEAD"])
        config["__worfklow_commit__"] = out.decode().strip()
except Exception as e:
    print(e)
    raise

##############################
## Config checks
##############################
check_blobdir_keys()

##################################################
## Uri parsing functions
##################################################
include: "core/uri.smk"


##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {
        'multiqc': [str(__REPORTS__ / "qc/multiqc.html")],
        'genecovr': all_genecovr_input(wildcards),
    }
    d['config'] = "config/assemblyeval.config.yaml"
    return d

##############################
# genecovr
##############################
include: "input/genecovr.smk"

##############################
# btk
##############################
include: "input/btk.smk"

##############################
# quast
##############################
include: "input/quast.smk"

##############################
# jellyfish
##############################
include: "input/jellyfish.smk"

##############################
# busco
##############################
include: "input/busco.smk"

##############################
# multiqc
##############################
include: "input/multiqc.smk"
