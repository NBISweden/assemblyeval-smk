import sys
import re
import os
import yaml
import itertools
import urllib
import pandas as pd
import numpy as np
import subprocess as sp

# Determine wrapper prefix since we mix local wrappers with wrappers
# from snakemake-wrappers
SMK_WRAPPER_VERSION = "v3.3.3"
SMK_WRAPPER_PREFIX_RAW = "https://github.com/snakemake/snakemake-wrappers/raw"
SMK_WRAPPER_PREFIX = f"{SMK_WRAPPER_PREFIX_RAW}/{SMK_WRAPPER_VERSION}"
try:
    WRAPPER_PREFIX = workflow.wrapper_prefix.rstrip("/")
except:
    WRAPPER_PREFIX = workflow.workflow_settings.wrapper_prefix.rstrip("/")

# Point wrapper prefix to github repo
if WRAPPER_PREFIX == SMK_WRAPPER_PREFIX_RAW:
    # Change main to version number once we start sem-versioning
    WRAPPER_PREFIX = "https://raw.githubusercontent.com/NBISweden/assemblyeval-smk/main/workflow/wrappers"


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##############################
# Core configuration
##############################
include: "core/config.smk"


##### load config and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

# If these are handled by config, try to get rid of them
# assemblies = Assemblies(config["assemblies"])
# transcripts = Transcripts(config["transcripts"])
# reads = Reads(config["reads"])

## Store some workflow metadata
config["__workflow_basedir__"] = workflow.basedir
config["__workflow_workdir__"] = os.getcwd()
config["__worfklow_commit__"] = None
config["__workflow_commit_link__"] = None

try:
    with cd(workflow.basedir, logger):
        commit = sp.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
        commit_short = (
            sp.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()
        )
        config["__workflow_commit__"] = commit_short
        config[
            "__workflow_commit_link__"
        ] = f"https://github.com/NBISweden/assemblyeval-smk/commit/{commit}"
except Exception as e:
    print(e)
    raise

# Wrap config dictionary
cfg = Config(config)

##############################
## Paths
##############################
__EXTERNAL__ = Path(cfg.fs.external)
__INTERIM__ = Path(cfg.fs.interim)
__METADATA__ = Path(cfg.fs.metadata)
__RAW__ = Path(cfg.fs.raw)
__REPORTS__ = Path(cfg.fs.reports)
__RESOURCES__ = Path(cfg.fs.resources)
__RESULTS__ = Path(cfg.fs.results)


##############################
## Wildcard constraints
##############################
wildcard_constraints:
    external=str(__EXTERNAL__),
    interim=str(__INTERIM__),
    metadata=str(__METADATA__),
    raw=str(__RAW__),
    reports=str(__REPORTS__),
    resources=str(__RESOURCES__),
    results=str(__RESULTS__),


## Assemblies etc
wildcard_constraints:
    analysis="|".join(cfg.analysisnames),
    assembly="|".join(cfg._assemblies.ids),
    blobdir="[^/]+",
    length="[0-9]+",
    partition="[0-9]+",
    region="[0-9]+",
    sep="(|/)",


## File extensions
wildcard_constraints:
    fa="(.fa|.fasta)",
    gz="(|.gz)",


##################################################
# Input collection functions
##################################################
def all(wildcards):
    d = {
        "multiqc": [str(__REPORTS__ / "multiqc.html")],
        "genecovr": all_genecovr(wildcards),
    }
    d["config"] = "config/assemblyeval.config.yaml"
    return d


##############################
# genecovr
##############################
include: "input/genecovr.smk"


##############################
# btk
##############################
# include: "input/btk.smk"


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
