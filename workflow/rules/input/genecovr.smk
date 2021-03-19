__GENECOVR_MINMATCH__ = [0.75, 0.85, 0.9, 0.95]
__GENECOVR_NCONTIGS__ = ["bar"]
__GENECOVR_MATCH_INDEL__ = ["violin", "boxplot", "boxplot.log10"]
__GENECOVR_FN__ = ["width_violin.pdf", "qnuminsert.pdf"]
__GENECOVR_DEPTH_BREADTH__ = ["coverage", "jitter", "hist", "seqlengths"]
__GENECOVR_CSV_GZ__ = ["gene_body_coverage.csv.gz", "psldata.csv.gz", "gbc_summary.csv.gz", "ncontigs_per_transcripts.csv.gz"]


def get_genecovr_input(wildcards):
    for analysis, toolconf, assembly_ids, read_ids, transcript_ids in iter_analyses("genecovr"):
        if analysis != wildcards.analysis:
            continue
        pfx = __RESULTS__ / f"genecovr/{analysis}/"
        outputs = genecovr_output(pfx)
        d = genecovr_make_csv_inputfile_input(assembly_ids, transcript_ids)
        csvfile = str(pfx / "genecovr.csv")
        if not os.path.exists(csvfile):
            if not os.path.exists(pfx):
                os.makedirs(pfx)
            df = pd.DataFrame(d)
            validate(df, schema=os.path.join(BASEDIR, "../schemas/genecovr_csv.schema.yaml"))
            df.to_csv(csvfile, index=False, header=False)
        break
    retval = d["psl"] + d["assembly"] + d["trxset"]
    return {'csv': csvfile, 'files': retval}


def genecovr_output(pfx=None):
    retval = []
    retval += report(expand(pfx / "gene_body_coverage.minmatch.{mm}.pdf", mm=__GENECOVR_MINMATCH__),
                     caption="../report/genecovr_gbc.rst", category="Gene body coverages")
    retval += report(expand(pfx / "ncontigs_per_transcripts.{type}.mm0.75.pdf", type=__GENECOVR_NCONTIGS__),
                     caption="../report/genecovr_ncontigs.rst", category="Number of contigs per transcript")
    retval += report(expand(pfx / "depth_breadth_{type}.mm0.75.pdf", type=__GENECOVR_DEPTH_BREADTH__),
                     caption="../report/genecovr_depth_breadth.rst", category="Depth and breadth of coverage")
    retval += report(expand(pfx / "match_indel.{type}.pdf", type=__GENECOVR_MATCH_INDEL__),
                     caption="../report/genecovr_match_indel.rst", category="Match and indel distributions")
    retval += report(expand(pfx / "{fn}", fn=__GENECOVR_FN__),
                     caption="../report/genecovr_match_indel.rst", category="Match and indel distributions")
    retval += report(expand(pfx / "{fn}", fn=__GENECOVR_CSV_GZ__),
                     caption="../report/genecovr_data.rst", category="Data files")
    return retval


def all_genecovr(wildcards):
    val = []
    for analysis, toolconf, assembly_ids, read_ids, transcript_ids in iter_analyses("genecovr"):
        pfx = __RESULTS__ / f"genecovr/{analysis}/"
        outputs = genecovr_output(pfx)
        val.extend(outputs)
    return val


def genecovr_make_csv_inputfile_input(assembly_ids, transcript_ids):
    assembly_fasta = assemblies.loc[assembly_ids]["fasta"].tolist()
    trx_fasta = transcripts.loc[transcript_ids, "fasta"].tolist()
    retval = {
        'dataset': [f"{a}/{b}" for a, b in itertools.product(assembly_ids, transcript_ids)],
        'psl': [str(__INTERIM__ / f"gmap/map/{a}-{b}.psl") for a, b in itertools.product(assembly_ids, transcript_ids)],
        'assembly': [f"{a}.fai" for a, b in itertools.product(assembly_fasta, trx_fasta)],
        'trxset': [f"{b}" for a, b in itertools.product(assembly_fasta, trx_fasta)]
    }
    return retval
