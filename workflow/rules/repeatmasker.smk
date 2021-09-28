rule repeatmasker_make_chunked_input:
    """Chunkify assembly input"""
    output:
        temp(
            expand(
                "{{interim}}/repeatmasker/{{assembly}}/fasta/{partition}.fasta",
                partition=range(config.get("repeatmasker", {}).get("npartitions", 1)),
            )
        ),
    input:
        seq=lambda wildcards: cfg.get_assembly(wildcards.assembly),
        faidx=lambda wildcards: cfg.get_assembly(wildcards.assembly, fai=True),
    conda:
        "../envs/pybedtools.yaml"
    resources:
        runtime=cfg.ruleconf("repeatmasker_make_chunked_input").xruntime,
    threads: cfg.ruleconf("repeatmasker_make_chunked_input").threads
    log:
        "logs/{interim}/repeatmasker/{assembly}/fasta/repeatmasker_make_chunked_input.log",
    script:
        "../scripts/assemblyeval_pybedtools_make_chunks.py"


# RepeatMasker -a -pa 20 -engine ncbi -gff -lib <library> <fasta_chunk>
rule repeatmasker_chunk:
    """Run RepeatMasker on chunk"""
    output:
        align=temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.align"),
        tbl=temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.tbl"),
        cat=temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.cat"),
        masked=temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.masked"),
        out=temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.out"),
    input:
        "{interim}/repeatmasker/{assembly}/fasta/{partition}.fasta",
    params:
        options=cfg.ruleconf("repeatmasker_chunk").options,
    envmodules:
        *cfg.ruleconf("repeatmasker_chunk").envmodules,
    resources:
        runtime=cfg.ruleconf("repeatmasker_chunk").xruntime,
        mem_mb=cfg.ruleconf("repeatmasker_chunk").xmem,
    log:
        "logs/{interim}/repeatmasker/{assembly}/{partition}.log",
    threads: cfg.ruleconf("repeatmasker_chunk").threads
    wrapper:
        f"{WRAPPER_PREFIX}/bio/repeatmasker"


rule repeatmasker_all:
    output:
        "results/repeatmasker/{assembly}/repeatmasker.ok",
    input:
        expand(
            __INTERIM__ / "repeatmasker/{{assembly}}/{partition}.fasta.align",
            partition=range(100),
        ),
    log:
        "logs/results/repeatmasker/{assembly}/repeatmasker.log",
    shell:
        "touch {output}"
