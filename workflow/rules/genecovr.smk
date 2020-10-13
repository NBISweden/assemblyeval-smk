wildcard_constraints:
    genecovr_results = str(__RESULTS__ / "genecovr")

__GENECOVR_MINMATCH__ = [0.75, 0.85, 0.9, 0.95]
__GENECOVR_NCONTIGS__ = ["bar"]
__GENECOVR_MATCH_INDEL__ = ["violin", "boxplot", "boxplot.log10"]
__GENECOVR_FN__ = ["width_violin.pdf", "qnuminsert.pdf"]
__GENECOVR_DEPTH_BREADTH__ = ["coverage", "dataset_jitter", "dataset_hist"]
__GENECOVR_CSV_GZ__ = ["gene_body_coverage.csv.gz", "psldata.csv.gz", "gbc_summary.csv.gz", "ncontigs_per_transcripts.csv.gz"]



rule genecovr_all:
    input: get_genecovr_all


def _genecovr_make_csv_inputfile_input(wildcards):
    assembly_keys = config["genecovr"][wildcards.dataset]["assemblies"]
    species = [k.split("_")[0] for k in assembly_keys]
    version = [k.split("_")[1] for k in assembly_keys]
    trx_keys = config["genecovr"][wildcards.dataset]["transcripts"]
    assembly_fasta = assemblies.loc[(species, version), "fasta"].tolist()
    trx_fasta = transcripts.loc[trx_keys, "fasta"].tolist()
    retval = {
        'psl': [str(__INTERIM__ / f"gmap/map/{a}-{b}.psl") for a, b in itertools.product(assembly_keys, trx_keys)],
        'assembly': [f"{a}" for a, b in itertools.product(assembly_fasta, trx_fasta)],
        'trxset': [f"{b}" for a, b in itertools.product(assembly_fasta, trx_fasta)]
    }
    return retval

rule genecovr_make_csv_inputfile:
    """Generate csv input file from dataset key"""
    wildcard_constraints:
        outprefix = "({})".format("|".join([config["genecovr"][x]["outprefix"] for x in config["genecovr"].keys() if x.startswith("dataset")]))
    output:
        csv = "{outprefix}.{dataset}.csv"
    input:
        unpack(_genecovr_make_csv_inputfile_input)
    run:
        assembly_keys = config["genecovr"][wildcards.dataset]["assemblies"]
        trx_keys = config["genecovr"][wildcards.dataset]["transcripts"]
        d = dict(input)
        df = pd.concat([
            pd.Series([f"{a}/{b}" for a, b in itertools.product(assembly_keys, trx_keys)]),
            pd.Series(list(d["psl"])), pd.Series(list(d["assembly"])),
            pd.Series(list(d["trxset"]))], axis=1)
        df.to_csv(output.csv, index=False, header=False)


def get_genecovr_input(wildcards):
    retval = []
    csvfile = config["genecovr"][wildcards.dataset]['csvfile']
    if csvfile is not None:
        data = pd.read_csv(csvfile, header=None).set_index(0, drop=False)
        data.columns = ["dataset", "psl", "assembly", "trxset"]
        for k in ["psl", "assembly", "trxset"]:
            retval += [x for x in data[k].tolist() if x is not None]
    else:
        csvfile = f"{config['genecovr'][wildcards.dataset]['outprefix']}.{wildcards.dataset}.csv"
    return {'csv': csvfile, 'files': retval}


def genecovr_output():
    retval = []
    retval += report(expand("{{genecovr_results}}/{{dataset}}/gene_body_coverage.minmatch.{mm}.pdf", mm=__GENECOVR_MINMATCH__),
                     caption="../report/genecovr.rst", category="Gene body coverages")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/ncontigs_per_transcripts.{type}.mm0.75.pdf", type=__GENECOVR_NCONTIGS__),
                     caption="../report/genecovr.rst", category="n contigs per transcripts")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/depth_breadth_{type}.mm0.75.pdf", type=__GENECOVR_DEPTH_BREADTH__),
                     caption="../report/genecovr.rst", category="Depth and breadth of coverage")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/match_indel.{type}.pdf", type=__GENECOVR_MATCH_INDEL__),
                     caption="../report/genecovr.rst", category="Match and indel distributions")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/{fn}", fn=__GENECOVR_FN__),
                     caption="../report/genecovr.rst", category="Match and indel distributions")
    retval += report(expand("{{genecovr_results}}/{{dataset}}/{fn}", fn=__GENECOVR_CSV_GZ__),
                     caption="../report/genecovr.rst", category="Data files")
    return retval


rule genecovr:
    """Calculate gene body coverage"""
    output:
        genecovr_output()
    input:
        unpack(get_genecovr_input)
    params:
        options = config["genecovr"]["options"],
        exe = config["genecovr"]["exe"]
    resources:
        runtime = lambda wildcards, attempt: attempt * config["genecovr"]["runtime"]
    conda:
        "../envs/genecovr.yaml"
    threads:
        config["genecovr"]["threads"]
    shell:
        "{params.exe} -p {threads} {params.options} -d results/genecovr/{wildcards.dataset} {input.csv}"


localrules: genecovr_make_csv_inputfile
