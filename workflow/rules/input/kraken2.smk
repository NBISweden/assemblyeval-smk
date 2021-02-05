def all_kraken2_input(wildcards):
    if "kraken2" not in config.keys():
        logger.info("kraken2 config section missing; no kraken2 analyses will be run")
        return []
    pfx = str(__RESULTS__ / "kraken2/{assembly}.{db}.{window_size}.report.txt")
    ids = make_assembly_ids(config["kraken2"].get("ids", []))
    val = expand(pfx, assembly=ids, window_size=config["kraken2"]["window_size"],
                 db=os.path.basename(config["kraken2"]["db"]))
    return val


def kraken2_gather_results_input(wildcards):
    fmt = __INTERIM__ / f"kraken2/{wildcards.assembly}/{wildcards.db}.{wildcards.length}.{{partition}}.{{suffix}}"
    d = {
        'output': expand(fmt, partition=range(0, get_workflow_params("kraken2", "npartitions")), suffix="output.txt.gz"),
        'unclassified': expand(fmt, partition=range(0, get_workflow_params("kraken2", "npartitions")), suffix="unclassified.fasta.gz")
    }
    return d


def kraken2_gather_reports_input(wildcards):
    fmt = __INTERIM__ / f"kraken2/{wildcards.assembly}/{wildcards.db}.{wildcards.length}.{{partition}}.report.txt"
    return expand(fmt, partition=range(0, get_workflow_params("kraken2", "npartitions")))
