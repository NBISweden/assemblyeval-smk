def all_kmer_comparison(wildcards):
    pfx = str(__RESULTS__ / "jellyfish/{assembly}.{analysis}.{kmer}_jf.png")
    val = []
    for analysis in config.keys():
        if not analysis.startswith("analysis"):
            continue
        jellyfish_plot = False
        for tool in config[analysis].get("tools", []):
            if tool["name"] == "jellyfish_plot":
                jellyfish_plot = True
                continue
        if jellyfish_plot:
            if "jellyfish" not in config.keys():
                logger.info("jellyfish config section missing but jellyfish_plot tool listed; jellyfish config to set kmer values")
                return []
            assembly = config[analysis]["ids"]
            val.extend(expand(pfx, analysis=[re.sub(r"^analysis/", "", analysis)],
                              assembly=assembly, kmer=config["jellyfish"]["kmer"]))
    return val


def all_jellyfish(wildcards):
    if "jellyfish" not in config.keys():
        logger.info("jellyfish config section missing; no jellyfish analyses will be run")
        return []
    pfx = str(__RESULTS__ / "jellyfish/assembly/{assembly}.{kmer}_jf.hist")
    kmer = get_workflow_params("jellyfish", "kmer")
    ids = make_assembly_ids(config["jellyfish"].get("ids", []))
    return expand(pfx, assembly=ids, kmer=kmer)


def jellyfish_count_input(wildcards):
    if wildcards.prefix in assemblies.index:
        wildcards.assembly = wildcards.prefix
        return {'seq': get_assembly(wildcards)}
    analysis = f"analysis/{wildcards.prefix}"
    try:
        read_ids = config[analysis].get("reads", [])
    except KeyError as e:
        print(e)
        raise
    reads = get_reads(read_ids)
    return {'seq': reads}
