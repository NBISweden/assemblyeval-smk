rule all_busco:
    input: all_busco_input


rule busco_run:
    output:
        tsv = "{results}/busco/{assembly}/{mode}/run_{lineage}/full_table.tsv",
        missing = "{results}/busco/{assembly}/{mode}/run_{lineage}/missing_busco_list.tsv",
        summary = "{results}/busco/{assembly}/{mode}/run_{lineage}/short_summary.txt"
    input:
        get_assembly
    wildcard_constraints:
        mode = "(genome|transcriptome|proteins)"
    resources:
        runtime = lambda wildcards, attempt: resources("busco_run", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("busco_run", "mem_mb", attempt),
    params:
        options = get_params("busco_run", "options")
    threads:
        get_params("busco_run", "threads")
    log: "logs/{results}/busco/{assembly}/{mode}/run_{lineage}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/busco"
