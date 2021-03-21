rule all_quast:
    input: unpack(all_quast)


rule quast:
    """Run quast"""
    output:
        tsv = "{results}/quast/{analysis}/{assembly}/report.tsv",
        tsv_trans = "{results}/quast/{analysis}/{assembly}/transposed_report.tsv",
        txt = "{results}/quast/{analysis}/{assembly}/report.txt",
        txt_trans = "{results}/quast/{analysis}/{assembly}/transposed_report.txt"
    input:
        seq = get_assembly
    resources:
        runtime = lambda wildcards, attempt: resources("quast", "runtime", attempt, wildcards.analysis),
        mem_mb = lambda wildcards, attempt: resources("quast", "mem_mb", attempt, wildcards.analysis),
    threads:
        lambda wildcards, attempt: resources("quast", "threads", attempt, wildcards.analysis)
    log:
        "logs/{results}/quast/{analysis}/{assembly}.log"
    envmodules: *get_params("quast", "envmodules")
    wrapper:
        f"{WRAPPER_PREFIX}/bio/quast/quast"
