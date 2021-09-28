def get_genecovr_input(wildcards):
    k = f"analysis/{wildcards.analysis}"
    genecovr = cfg[k].tools["genecovr"]
    d = genecovr.csv_input()
    csvfile = str(genecovr.csvfile).format(**wildcards)
    if not os.path.exists(csvfile):
        if not os.path.exists(os.path.dirname(csvfile)):
            os.makedirs(os.path.dirname(csvfile))
        df = pd.DataFrame(d)
        genecovr_schema.validate(df)
        df.to_csv(csvfile, index=False, header=False)
    retval = d["psl"] + d["assembly"] + d["trxset"]
    return {"csv": csvfile, "files": retval}


def all_genecovr(wildcards):
    val = []
    for analysis in cfg.analyses_w_tool("genecovr"):
        val.extend(analysis.tools["genecovr"].input)
    return val
