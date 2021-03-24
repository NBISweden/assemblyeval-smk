import contextlib

##############################
# Get configuration options
##############################
def get_toolconf(tool, option, analysis="default"):
    """Retrieve tool configuration

    This function should only be called if there is a rule that can
    generate results, i.e. the values should be set.

    """
    akey = f"analysis/{analysis}"
    if analysis == "default":
        toolconf = config["tools"][tool]
    else:
        toolconf = config[akey]["tools"][tool]
    return toolconf.get(option)


def resources(rule, resource, attempt=1, analysis="default", wildcards=None, **kwargs):
    """Retrieve resources for rule multiplying the value by attempt"""
    akey = f"analysis/{analysis}"
    subconf = config.get(akey, config)
    if subconf.get('rules', {}).get(rule, {}).get(resource, None):
        val = subconf['rules'][rule][resource]
    else:
        val = config['resources.default'][resource]

    return attempt * val


def get_params(rule, resource, wildcards=None, **kwargs):
    """Retrieve rule parameters from top-level configuration"""
    val = None
    if "rules" in config.keys():
        val = config["rules"].get(rule, {}).get(resource, None)
    if val is not None:
        return val
    val = config['resources.default'][resource]
    return val


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
    ids = ["default"] + [k.replace("analysis/", "") for k in config.keys() if k.startswith("analysis/")]
    return ids


def make_read_ids(ids=None):
    if len(ids) == 0 or ids is None:
        return reads.index.to_list()
    try:
        assert set(ids) <= set(reads.index.to_list())
    except AssertionError as e:
        logger.error(f"undefined read identifiers in ids: \'{', '.join(ids)}\'")
        raise
    return ids


def make_transcript_ids(ids=None):
    if len(ids) == 0 or ids is None:
        return transcripts.index.to_list()
    try:
        assert set(ids) <= set(transcripts.index.to_list())
    except AssertionError as e:
        logger.error(f"undefined transcript identifiers in ids: \'{', '.join(ids)}\'")
        raise
    return ids


def iter_analyses(tool):
    """Iterate analyses and retrieve tool configuration, assembly, transcript, and read ids"""
    for k, analysis in config.items():
        if not k.startswith("analysis/") and not k == "tools":
            continue
        # Top-level analysis definition
        if tool in config[k].keys():
            analysis = config
            k = "default"
        if tool not in analysis.get("tools", {}).keys():
            continue
        toolconf = analysis["tools"][tool]
        assembly_ids = make_assembly_ids(analysis.get("assembly_ids", []))
        read_ids = make_read_ids(analysis.get("read_ids", []))
        transcript_ids = make_transcript_ids(analysis.get("transcript_ids", []))
        yield k.replace("analysis/", ""), toolconf, assembly_ids, read_ids, transcript_ids



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
