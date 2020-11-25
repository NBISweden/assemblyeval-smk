def all_busco_input(wildcards):
    if "busco" not in config.keys():
        logger.info("busco config section missing; no busco analyses will be run")
        return []
    pfx = str(__RESULTS__ / "busco/{assembly}/{mode}/run_{lineage}/full_table.tsv")
    ids = make_assembly_ids(config["quast"].get("ids", []))
    val = expand(pfx, assembly=ids, mode=config["busco"].get("mode", ["genome"]),
                 lineage=config["busco"]["lineage"])
    return val
