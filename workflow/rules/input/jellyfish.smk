def all_kmer_comparison(wildcards):
    pfx = str(__RESULTS__ / "jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.png")
    val = []
    for analysis in config.keys():
        if not analysis.startswith("analysis"):
            continue
        jellyfish_plot = False
        for tool in config[analysis].get("tools", []):
            if tool["name"] == "jellyfish_kmer_pairs_plot":
                jellyfish_plot = True
                continue
        if jellyfish_plot:
            if "jellyfish" not in config.keys():
                logger.info("jellyfish config section missing but jellyfish_kmer_pairs_plot tool listed; jellyfish config to set kmer values")
                return []
            assembly = config[analysis]["ids"]
            val.extend(expand(pfx, analysis=[re.sub(r"^analysis/", "", analysis)],
                              assembly=assembly, kmer=config["jellyfish"]["kmer"]))
    return val


def all_jellyfish(wildcards):
    if "jellyfish" not in config.keys():
        logger.info("jellyfish config section missing; no jellyfish analyses will be run")
        return []
    pfx = str(__RESULTS__ / "jellyfish/{assembly}/merged.{kmer}_jf.hist")
    kmer = get_workflow_params("jellyfish", "kmer")
    ids = make_assembly_ids(config["jellyfish"].get("ids", []))
    return expand(pfx, assembly=ids, kmer=kmer)


def jellyfish_count_input(wildcards):
    if wildcards.dataset in assemblies.index:
        return {'seq': str(__INTERIM__ / "jellyfish/{dataset}/{prefix}.fasta".format(**wildcards))}
    read = get_reads(name=wildcards.prefix)
    return {'seq': read}


def jellyfish_merge_input(wildcards):
    if wildcards.dataset in assemblies.index:
        npartitions = config.get("jellyfish", {}).get("npartitions", 1)
        return {'seq': expand(str(__INTERIM__ / "jellyfish/{dataset}/{{partition}}.{kmer}mer_counts.jf".format(**wildcards)),
                              partition=range(npartitions))}
    analysis = f"analysis/{wildcards.dataset}"
    try:
        read_ids = config[analysis].get("reads", [])
    except KeyError as e:
        print(f"No such key '{e}' in config; using all reads")
        raise
    reads = get_reads(read_ids)
    return {'seq': expand(str(__INTERIM__ / "jellyfish/{dataset}/{{prefix}}.{kmer}mer_counts.jf".format(**wildcards)),
                          prefix=[os.path.basename(x) for x in reads])}
