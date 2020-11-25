def all_jellyfish(wildcards):
    if "jellyfish" not in config.keys():
        logger.info("jellyfish config section missing; no jellyfish analyses will be run")
        return []
    pfx = str(__RESULTS__ / "jellyfish/{assembly}.{kmer}_jf.hist")
    kmer = get_workflow_params("jellyfish", "kmer")
    ids = make_assembly_ids(config["jellyfish"].get("ids", []))
    return expand(pfx, assembly=ids, kmer=kmer)
