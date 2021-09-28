def all_quast(wildcards):

    val = {"tsv": [], "other": []}
    for analysis in cfg.analyses_w_tool("quast"):
        files = analysis.tools["quast"].input
        val["tsv"].extend(filter(lambda fn: re.search(r".*/report.tsv$", fn), files))
        val["other"].extend(
            filter(lambda fn: not re.search(r".*/report.tsv$", fn), files)
        )
    return val
