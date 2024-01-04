def all_multiqc(wildcards):
    val = {
        "quast": all_quast(wildcards).get("tsv", []),
        "jellyfish": all_jellyfish(wildcards),
        "busco": all_busco_input(wildcards),
        "kraken": all_kraken2_input(wildcards),
    }
    return val
