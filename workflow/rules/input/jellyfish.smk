def all_jellyfish(wildcards):
    val = []
    for analysis, toolconf, assembly_ids, _, _ in iter_analyses("jellyfish"):
        pfx = __RESULTS__ / f"jellyfish/{analysis}/{{assembly}}/merged.{{kmer}}_jf.hist"
        val.extend(expand(pfx, assembly=assembly_ids, kmer=toolconf["kmer"]))
        if toolconf["count_pairs"]:
            pfx = __RESULTS__ / f"jellyfish/{analysis}/kmer_comparison/{{assembly}}.{{kmer}}_jf.png"
            val.extend(expand(pfx, assembly=assembly_ids, kmer=toolconf["kmer"]))
    return val


def jellyfish_count_input(wildcards):
    """Get input file for jellyfish count. Return either an assembly file or read files"""
    if wildcards.dataset in assemblies.index:
        return {'seq': str(__INTERIM__ / "jellyfish/{analysis}/{dataset}/{prefix}.fasta".format(**wildcards))}
    read = get_reads(name=wildcards.prefix)
    return {'seq': read}


def jellyfish_merge_input(wildcards):
    if wildcards.dataset in assemblies.index:
        npartitions = get_toolconf("jellyfish", "npartitions", wildcards.analysis)
        return {'jf': expand(str(__INTERIM__ / "jellyfish/{analysis}/{dataset}/{{partition}}.{kmer}mer_counts.jf".format(**wildcards)),
                             partition=range(npartitions))}
    analysis = f"analysis/{wildcards.dataset}"
    read_ids = []
    try:
        read_ids = config[analysis].get("reads", [])
    except KeyError as e:
        print(f"No such key '{e}' in config; using all reads")
    reads = get_reads(read_ids)
    return {'jf': expand(str(__INTERIM__ / "jellyfish/{analysis}/{dataset}/{{prefix}}.{kmer}mer_counts.jf".format(**wildcards)),
                         prefix=[os.path.basename(x) for x in reads])}
