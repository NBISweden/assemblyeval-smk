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
    val = config['resources'][rule].get(resource, None)
    if val is not None:
        return val
    val = config['resources.default'][resource]
    return val


def check_blobdir_keys():
    """Check that blobdir key is formatted correctly and that species and version exist"""
    blobdir_assemblies = []
    for blobdir in config["btk"].keys():
        if not blobdir.startswith("blobdir_"):
            continue
        species, version = re.sub("blobdir_", "", blobdir).split("_")
        if (species, version) not in assemblies.index:
            logger.error(f"error in blobdir configuration: {species}_{version} not in {config['assemblies']}")
            sys.exit(1)


##################################################
## Getters for assembly and transcriptome
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


# context manager for cd
@contextlib.contextmanager
def cd(path, logger):
    CWD = os.getcwd()
    logger.info("Changing directory from {} to {}".format(CWD, path))

    os.chdir(path)
    try:
        yield
    except:
        logger.warning("Exception caught: ".format(sys.exc_info()[0]))
    finally:
        logger.info("Changing directory back to {}".format(CWD))
        os.chdir(CWD)
