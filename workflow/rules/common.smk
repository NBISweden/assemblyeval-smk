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
SMK_WRAPPER_PREFIX = "https://github.com/snakemake/snakemake-wrappers/raw/"
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
## Uri parsing functions
##################################################
def get_uri_scheme(uri):
    return urllib.parse.urlparse(uri).scheme

def get_uri_netloc(uri):
    return urllib.parse.urlparse(uri).netloc

def parse_uri(uri):
    """Parse uri and return snakemake target"""
    scheme = get_uri_scheme(uri)
    uri = re.sub(f"{scheme}://", "", uri)
    if not scheme in ['', 'rsync', 'file', 'sftp']:
        logger.error(f"scheme '{scheme}' not allowed: use one of '', 'file', 'rsync' or 'sftp'")
        sys.exit(1)
    if scheme in ['', 'file'] and not uri.startswith("/"):
        uri = os.path.normpath(os.path.abspath(uri))
    if scheme == 'sftp':
        try:
            from snakemake.remote.SFTP import RemoteProvider
            SFTP = RemoteProvider()
            uri = SFTP.remote(uri)
        except WorkflowError as e:
            logger.error(e)
    return uri

def assemblyeval_get_external_input(uri):
    if get_uri_scheme(uri) == "sftp":
        return parse_uri(uri)
    netloc = get_uri_netloc(uri)
    uri = parse_uri(uri)
    if re.search("[:@]", netloc):
        return []
    return uri

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


def genecovr_make_csv_inputfile_input(wildcards):
    assembly_keys = config["genecovr"][wildcards.dataset]["assemblies"]
    species = [k.split("_")[0] for k in assembly_keys]
    version = [k.split("_")[1] for k in assembly_keys]
    trx_keys = config["genecovr"][wildcards.dataset]["transcripts"]
    assembly_fasta = assemblies.loc[(species, version), "fasta"].tolist()
    trx_fasta = transcripts.loc[trx_keys, "fasta"].tolist()
    retval = {
        'psl': [str(__INTERIM__ / f"gmap/map/{a}-{b}.psl") for a, b in itertools.product(assembly_keys, trx_keys)],
        'assembly': [f"{a}.fai" for a, b in itertools.product(assembly_fasta, trx_fasta)],
        'trxset': [f"{b}" for a, b in itertools.product(assembly_fasta, trx_fasta)]
    }
    return retval

def genecovr_make_csv_inputfile_dataframe(d, wildcards):
    assembly_keys = config["genecovr"][wildcards.dataset]["assemblies"]
    trx_keys = config["genecovr"][wildcards.dataset]["transcripts"]
    df = pd.concat([
        pd.Series([f"{a}/{b}" for a, b in itertools.product(assembly_keys, trx_keys)]),
        pd.Series(list(d["psl"])), pd.Series(list(d["assembly"])),
        pd.Series(list(d["trxset"]))], axis=1)
    return df


def get_genecovr_input(wildcards):
    retval = []
    csvfile = config["genecovr"][wildcards.dataset]['csvfile']
    if csvfile is not None:
        data = pd.read_csv(csvfile, header=None).set_index(0, drop=False)
        data.columns = ["dataset", "psl", "assembly", "trxset"]
        try:
            assert data["assembly"].str.replace(r".fai", "").isin(assemblies["fasta"]).all()
        except AssertionError as e:
            logger.error(e)
            logger.error("some values in 'assembly' column not present in assemblies input file")
        try:
            assert data["trxset"].str.replace(r".fai", "").isin(transcripts["fasta"]).all()
        except AssertionError as e:
            logger.error(e)
            logger.error("some values in 'trxset' column not present in transcripts input file")
        validate(data, schema=os.path.join(BASEDIR, "../schemas/genecovr_csv.schema.yaml"))
        for k in ["psl", "assembly", "trxset"]:
            retval += [x for x in data[k].tolist() if x is not None]
    else:
        csvfile = f"{config['genecovr'][wildcards.dataset]['outprefix']}.{wildcards.dataset}.csv"
    return {'csv': csvfile, 'files': retval}


def genecovr_output():
    retval = []
    retval += report(expand("{{genecovr_results}}/{{dataset}}/gene_body_coverage.minmatch.{mm}.pdf", mm=__GENECOVR_MINMATCH__),
                     caption="../report/genecovr_gbc.rst", category="Gene body coverages")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/ncontigs_per_transcripts.{type}.mm0.75.pdf", type=__GENECOVR_NCONTIGS__),
                     caption="../report/genecovr_ncontigs.rst", category="Number of contigs per transcript")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/depth_breadth_{type}.mm0.75.pdf", type=__GENECOVR_DEPTH_BREADTH__),
                     caption="../report/genecovr_depth_breadth.rst", category="Depth and breadth of coverage")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/match_indel.{type}.pdf", type=__GENECOVR_MATCH_INDEL__),
                     caption="../report/genecovr_match_indel.rst", category="Match and indel distributions")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/{fn}", fn=__GENECOVR_FN__),
                     caption="../report/genecovr_match_indel.rst", category="Match and indel distributions")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/{fn}", fn=__GENECOVR_CSV_GZ__),
                     caption="../report/genecovr_data.rst", category="Data files")
    return retval


def _btk_link_fasta_input(wildcards):
    return [_btk_link_fasta_input_paths(wildcards)['path']]


def _btk_link_fasta_input_paths(wildcards):
    ret = {}
    fn = Path(wildcards.interim) / "btk" / wildcards.blobdir / f"{wildcards.prefix}.fasta.gz"
    for x in config["btk"][wildcards.blobdir]["fasta"]:
        if str(fn.name) == str(Path(x).name):
            ret['path'] = str(Path(x))
            ret['abspath'] = str(Path(x).absolute())
            break
    return ret
