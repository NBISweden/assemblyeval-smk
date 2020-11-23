rule all_quast:
    input: unpack(all_quast)


rule quast:
    """Run quast"""
    output:
        tsv = "{results}/quast/{assembly}/report.tsv",
        tsv_trans = "{results}/quast/{assembly}/transposed_report.tsv",
        txt = "{results}/quast/{assembly}/report.txt",
        txt_trans = "{results}/quast/{assembly}/transposed_report.txt"
    input:
        seq = get_assembly
    resources:
        runtime = lambda wildcards, attempt: resources("quast", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("quast", "mem_mb", attempt),
    threads:
        lambda wildcards, attempt: resources("quast", "threads", attempt)
    log:
        "logs/{results}/quast/{assembly}.log"
    wrapper:
        f"{WRAPPER_PREFIX}/bio/quast/quast"
