def all_busco_input(wildcards):
    if "busco" not in config.keys():
        logger.info("busco config section missing; no busco analyses will be run")
        return []
    pfx = str(__RESULTS__ / "busco/{assembly}/{mode}/run_{lineage}/short_summary_{assembly}.txt")
    ids = make_assembly_ids(config["busco"].get("ids", []))
    val = expand(pfx, assembly=ids, mode=config["busco"].get("mode", ["genome"]),
                 lineage=config["busco"]["lineage"])
    return val
