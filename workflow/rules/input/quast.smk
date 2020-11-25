def all_quast(wildcards):
    if "quast" not in config.keys():
        logger.info("quast config section missing; no quast analyses will be run")
        return {}
    pfx = str(__RESULTS__ / "quast/{assembly}/{rpt}")
    rpt = ["report.tsv", "transposed_report.tsv", "report.txt",
           "transposed_report.txt"]
    ids = make_assembly_ids(config["quast"].get("ids", []))
    val = {
        'tsv': expand(pfx, assembly=ids, rpt="report.tsv"),
        'tsv.trans': expand(pfx, assembly=ids, rpt="transposed_report.tsv"),
        'txt': expand(pfx, assembly=ids, rpt="report.txt"),
        'txt.trans': expand(pfx, assembly=ids, rpt="report.txt")
    }
    return val
