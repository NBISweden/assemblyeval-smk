rule all_busco:
    input:
        all_busco_input,


rule busco_get_file_versions:
    """Retrieve busco file versions"""
    output:
        "busco_downloads/file_versions.tsv",
    envmodules:
        *cfg.ruleconf("busco_run").envmodules,
    log:
        "logs/busco_get_file_versions/busco_downloads/file_versions.tsv.log",
    threads: 1
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/busco/fv")


rule busco_download:
    """Download busco database"""
    output:
        "busco_downloads/lineages/{lineage}/dataset.cfg",
    input:
        "busco_downloads/file_versions.tsv",
    params:
        url=busco_lineage_input,
    log:
        "logs/busco_download/{lineage}.log",
    threads: 1
    shell:
        "wget {params.url} -P busco_downloads/lineages && tar -C busco_downloads/lineages -zxvf busco_downloads/lineages/$(basename {params.url});"
        "rm -f busco_downloads/lineages/$(basename {params.url})"


rule busco_run:
    output:
        tsv="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/full_table.tsv",
        missing="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/missing_busco_list.tsv",
        summary="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/short_summary.txt",
        multiqc_summary="{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}/short_summary_{assembly}.txt",
    input:
        assembly=lambda wildcards: cfg.get_assembly(wildcards.assembly),
        cfg="busco_downloads/lineages/{lineage}/dataset.cfg",
    wildcard_constraints:
        mode="(genome|transcriptome|proteins)",
    params:
        options=cfg.ruleconf("busco_run").options,
    envmodules:
        *cfg.ruleconf("busco_run").envmodules,
    threads: cfg.ruleconf("busco_run").xthreads
    log:
        "logs/{results}/busco/{analysis}/{assembly}/{mode}/run_{lineage}.log",
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/busco/run")


localrules:
    busco_get_file_versions,
