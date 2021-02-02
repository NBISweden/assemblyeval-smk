rule repeatmasker_make_chunked_input:
    """Chunkify assembly input"""
    output:
        temp(expand("{{interim}}/repeatmasker/{{assembly}}/fasta/{partition}.fasta", partition=range(config.get("repeatmasker", {}).get("npartitions", 1))))
    input:
        seq = get_assembly,
        faidx = get_assembly_index
    conda:
        "../envs/pybedtools.yaml"
    resources:
        runtime = lambda wildcards: resources("repeatmasker_make_chunked_input", "runtime")
    threads:
        get_params("repeatmasker_make_chunked_input", "threads")
    log: "logs/{interim}/repeatmasker/{assembly}/fasta/repeatmasker_make_chunked_input.log"
    script:
        "../scripts/assemblyeval_pybedtools_make_chunks.py"


# RepeatMasker -a -pa 20 -engine ncbi -gff -lib <library> <fasta_chunk>
rule repeatmasker_chunk:
    """Run RepeatMasker on chunk"""
    output:
        align = temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.align"),
        tbl = temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.tbl"),
        cat = temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.cat"),
        masked = temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.masked"),
        out = temp("{interim}/repeatmasker/{assembly}/{partition}.fasta.out")
    input:
        "{interim}/repeatmasker/{assembly}/fasta/{partition}.fasta"
    params:
        options = get_params("repeatmasker_chunk", "options")
    #envmodules: "\n".join([f"\"{x}\"" for x in get_params("repeatmasker_chunk", "envmodules")])
    envmodules:
        "bioinfo-tools",
        "RepeatMasker"
    resources:
        runtime = lambda wildcards, attempt: resources("repeatmasker_chunk", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("repeatmasker_chunk", "mem_mb", attempt),
    log: "logs/{interim}/repeatmasker/{assembly}/{partition}.log"
    threads:
        get_params("repeatmasker_chunk", "threads")
    wrapper:
        f"{WRAPPER_PREFIX}/bio/repeatmasker"


rule repeatmasker_all:
    output: "results/repeatmasker/{assembly}/repeatmasker.ok"
    input: expand(__INTERIM__ / "repeatmasker/{{assembly}}/{partition}.fasta.align", partition=range(100))
    shell:
        "touch {output}"
