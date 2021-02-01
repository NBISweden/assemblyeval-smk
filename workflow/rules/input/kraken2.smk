def all_kraken2_input(wildcards):
    if "kraken2" not in config.keys():
        logger.info("kraken2 config section missing; no kraken2 analyses will be run")
        return []
    pfx = str(__RESULTS__ / "kraken2/{assembly}.{db}.{window_size}.report.txt")
    ids = make_assembly_ids(config["kraken2"].get("ids", []))
    val = expand(pfx, assembly=ids, window_size=config["kraken2"]["window_size"],
                 db=os.path.basename(config["kraken2"]["db"]))
    return val
