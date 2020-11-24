def all_genecovr_input(wildcards):
    dataset = [x for x in config["genecovr"].keys() if x.startswith("dataset")]
    return expand(f"{str(__RESULTS__)}/genecovr/{{dataset}}/psldata.csv.gz", dataset=dataset)


def genecovr_make_csv_inputfile_input(wildcards):
    assembly_keys = config["genecovr"][wildcards.dataset]["assemblies"]
    trx_keys = config["genecovr"][wildcards.dataset]["transcripts"]
    assembly_fasta = assemblies.loc[assembly_keys]["fasta"].tolist()
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
