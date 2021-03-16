import contextlib


def resources(rule, resource, attempt=1, wildcards=None, **kwargs):
    """Retrieve resources for rule multiplying the value by attempt"""
    # TODO: Add prior checks to resource
    if config['resources'][rule].get(resource, None):
        val = config['resources'][rule][resource]
    else:
        val = config['resources.default'][resource]

    return attempt * val


def get_params(rule, resource, wildcards=None, **kwargs):
    """Retrieve rule parameters"""
    val = config["resources"][rule].get(resource, None)
    if val is not None:
        return val
    val = config['resources.default'][resource]
    return val


def get_workflow_params(section, option, **kwargs):
    """Retrieve workflow settings.

    This function works on parameters that cannot be set in the analysis sections"""
    if section not in config.keys():
        return
    val = config[section].get(option, None)
    return val


def check_blobdir_keys():
    """Check that blobdir key is formatted correctly and that species and version exist"""
    if "btk" not in config.keys():
        return
    blobdir_assemblies = []
    for blobdir in config["btk"].keys():
        if not blobdir.startswith("blobdir_"):
            continue
        blobid = re.sub("blobdir_", "", blobdir)
        if blobid not in assemblies.index:
            logger.error(f"error in blobdir configuration: key {blobid} not in {config['assemblies']}")
            sys.exit(1)


##################################################
## Getters for assembly and transcriptome
##################################################
def make_assembly_ids(ids=[]):
    """Make a complete list of assembly ids

    :param list ids: assembly identifier list

    :return list ids:  assembly identifier list
    """
    if len(ids) == 0:
        return assemblies.index.to_list()
    try:
        assert set(ids) <= set(assemblies.index.to_list())
    except AssertionError as e:
        logger.error(f"undefined assembly identifiers in ids: \'{', '.join(ids)}\'")
        raise
    return ids


def make_analysis_ids():
    """Make a complete list of analysis ids"""
    ids = [k.lstrip("analysis/") for k in config.keys() if k.startswith("analysis/")]
    return ids


def get_assembly(wildcards):
    """Retrieve the sequence file for a given assembly id"""
    # In principle could use ensembl wrapper if db not present:
    # file:// would be better
    return assemblies.loc[wildcards.assembly]["fasta"]


def get_assembly_index(wildcards):
    """Retrieve the sequence file for a given assembly id"""
    # In principle could use ensembl wrapper if db not present:
    # file:// would be better
    return assemblies.loc[wildcards.assembly]["fasta"] + ".fai"


def get_transcriptome(wildcards):
    """Retrieve the sequence file for a given transcriptome id"""
    value = transcripts.loc[wildcards.transcriptome]["fasta"]
    if isinstance(value, str):
        value = [value]
    return value


def get_reads(ids=None, name=None):
    """Retrieve the sequence files for a set of read ids"""
    allreads = reads["read1"].to_list() + reads["read2"].dropna().to_list()
    if name is not None:
        r = re.compile(f".+/{name}$")
        readlist = list(filter(r.search, allreads))
        return readlist
    if ids is None or len(ids) == 0:
        return allreads
    return reads.loc[ids]["read1"].to_list() + reads.loc[ids]["read2"].dropna().to_list()


# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except Exception as e:
        logger.warning(e)
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)
