def all_kraken2_input(wildcards):
    val = []
    for analysis, toolconf, assembly_ids, _, _ in iter_analyses("kraken2"):
        pfx = __RESULTS__ / "kraken2/{analysis}/{assembly}/{db}.{window_size}.report.txt"
        ids = make_assembly_ids(assembly_ids)
        tmp = expand(pfx, analysis=analysis, assembly=ids,
                     window_size=get_toolconf("kraken2",
                                              "window_size"),
                     db=os.path.basename(toolconf["db"]))
        val.extend(tmp)
    return val


def kraken2_gather_results_input(wildcards):
    fmt = __INTERIM__ / f"kraken2/{wildcards.analysis}/{wildcards.assembly}/{wildcards.db}.{wildcards.length}.{{partition}}.{{suffix}}"
    d = {
        'output': expand(fmt, partition=range(0, get_toolconf("kraken2", "npartitions", wildcards.analysis)), suffix="output.txt.gz"),
        'unclassified': expand(fmt, partition=range(0, get_toolconf("kraken2", "npartitions", wildcards.analysis)), suffix="unclassified.fasta.gz")
    }
    return d


def kraken2_gather_reports_input(wildcards):
    fmt = __INTERIM__ / f"kraken2/{wildcards.analysis}/{wildcards.assembly}/{wildcards.db}.{wildcards.length}.{{partition}}.report.txt"
    return expand(fmt, partition=range(0, get_toolconf("kraken2", "npartitions", wildcards.analysis)))
