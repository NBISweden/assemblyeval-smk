rule all_quast:
    input:
        unpack(all_quast),


rule quast:
    """Run quast"""
    output:
        tsv="{results}/quast/{analysis}/{assembly}/report.tsv",
        tsv_trans="{results}/quast/{analysis}/{assembly}/transposed_report.tsv",
        txt="{results}/quast/{analysis}/{assembly}/report.txt",
        txt_trans="{results}/quast/{analysis}/{assembly}/transposed_report.txt",
    input:
        seq=lambda wildcards: cfg.analysis(wildcards.analysis).get_assembly(
            wildcards.assembly
        ),
    threads: cfg.ruleconf("quast").xthreads
    log:
        "logs/{results}/quast/{analysis}/{assembly}.log",
    envmodules:
        *cfg.ruleconf("quast").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/quast/quast")
