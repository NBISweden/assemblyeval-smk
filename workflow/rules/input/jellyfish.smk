def all_jellyfish(wildcards):
    if "jellyfish" not in config.keys():
        logger.info("jellyfish config section missing; no jellyfish analyses will be run")
        return []
    pfx = str(__RESULTS__ / "jellyfish/{assembly}.{kmer}_jf.hist")
    kmer = get_workflow_params("jellyfish", "kmer")
    return expand(pfx, assembly=make_assembly_ids(), kmer=kmer)
