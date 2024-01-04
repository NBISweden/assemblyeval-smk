rule all_busco:
    input:
        all_busco_input,


rule busco_run:
    output:
        tsv="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/full_table.tsv",
        missing="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/missing_busco_list.tsv",
        summary="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/short_summary.txt",
        multiqc_summary="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/short_summary_{assembly}.txt",
    input:
        lambda wildcards: cfg.get_assembly(wildcards.assembly),
    wildcard_constraints:
        mode="(genome|transcriptome|proteins)",
    resources:
        runtime=cfg.ruleconf("busco_run").xruntime,
        mem_mb=cfg.ruleconf("busco_run").xmem,
    params:
        options=cfg.ruleconf("busco_run").options,
    envmodules:
        *cfg.ruleconf("busco_run").envmodules,
    threads: cfg.ruleconf("busco_run").xthreads
    log:
        "logs/{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}.log",
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/busco")
