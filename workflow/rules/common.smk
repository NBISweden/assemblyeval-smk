import sys
import re
import os
import yaml
import itertools
import urllib
import pandas as pd
import numpy as np
import subprocess as sp
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
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/NBISweden/assemblyeval-smk/main/workflow/wrappers"

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")


def _read(infile, index, schema, idcols=None):
    if infile is None:
        return None
    if os.path.splitext(infile)[1] == ".yaml":
        with open(infile) as fh:
            data = yaml.load(fh, yaml.Loader)
        assert(isinstance(data, list))
        df = pd.DataFrame(data)
    elif os.path.splitext(infile)[1] == ".tsv":
        df = pd.read_csv(infile, sep="\t")
    if "id" not in df.columns and index == ["id"]:
        logger.info(f"generating id column from {idcols}")
        df["id"] = "_".join(df[idcols])
        df["id"] = df[idcols].agg('_'.join, axis=1)
    df.set_index(index, drop=False, inplace=True)
    df = df.replace({np.nan: None})
    df.index.names = index
    validate(df, schema=schema)
    return df


assemblies = _read(config["assemblies"], ["id"],
                   "../schemas/assemblies.schema.yaml", idcols=["species", "version"])
datasources = _read(config["datasources"], ["data"],
                        "../schemas/datasources.schema.yaml")
transcripts = _read(config["transcripts"], ["id"],
                        "../schemas/transcripts.schema.yaml")
reads = _read(config["reads"], ["id"],
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
    analysis = "|".join(make_analysis_ids()),
    assembly = "|".join(make_assembly_ids()),
    blobdir = "[^/]+",
    length = "[0-9]+",
    partition = "[0-9]+",
    region = "[0-9]+",
    sep = "(|/)"

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

##################################################
## Uri parsing functions
##################################################
include: "core/uri.smk"

##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {
        'multiqc': [str(__REPORTS__ / "multiqc.html")],
        'genecovr': all_genecovr(wildcards)
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
# kraken2
##############################
include: "input/kraken2.smk"

##############################
# multiqc
##############################
include: "input/multiqc.smk"
