rule jellyfish_make_chunked_input:
    """Chunkify assembly input"""
    output:
        temp(expand("{{interim}}/jellyfish/{{assembly}}/{partition}.fasta", partition=range(config.get("jellyfish", {}).get("npartitions", 1))))
    input:
        seq = get_assembly,
        faidx = get_assembly_index
    conda:
        "../envs/pybedtools.yaml"
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_make_chunked_input", "runtime", attempt)
    threads:
        get_params("jellyfish_make_chunked_input", "threads")
    log: "logs/{interim}/jellyfish/{assembly}/fasta/jellyfish_make_chunked_input.log"
    script:
        "../scripts/assemblyeval_pybedtools_make_chunks.py"


rule jellyfish_count_chunk:
    """Count kmers in chunk"""
    output: jf = temp("{interim}/jellyfish/{dataset}/{prefix}.{kmer}mer_counts.jf")
    input: unpack(jellyfish_count_input)
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_count_chunk", "runtime", attempt)
    params:
        options = get_params("jellyfish_count_chunk", "options"),
        tmpdir = get_workflow_params("jellyfish", "tmpdir")
    threads: lambda wildcards, attempt: resources("jellyfish_count_chunk", "threads", attempt)
    log: "logs/{interim}/jellyfish/{dataset}/{prefix}.{kmer}mer_counts.log"
    envmodules: *get_params("jellyfish_count_chunk", "envmodules")
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/count"


rule jellyfish_merge:
    """Merge chunked kmer counts"""
    output: jf = "{results}/jellyfish/{dataset}/merged.{kmer}mer_counts.jf"
    input: unpack(jellyfish_merge_input)
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_merge", "runtime", attempt)
    params:
        options = get_params("jellyfish_merge", "options"),
        tmpdir = get_workflow_params("jellyfish", "tmpdir")
    log: "logs/{results}/jellyfish/{dataset}/merged.{kmer}mer_counts.jf.log"
    envmodules: *get_params("jellyfish_merge", "envmodules")
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/merge"



rule jellyfish_histo:
    output: hist = "{prefix}.{kmer}_jf.hist"
    input: counts = "{prefix}.{kmer}mer_counts.jf"
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_histo", "runtime", attempt)
    params:
        options = get_params("jellyfish_histo", "options"),
        tmpdir = get_workflow_params("jellyfish", "tmpdir")
    threads: get_params("jellyfish_histo", "threads")
    log: "logs/{prefix}.{kmer}_jf.hist.log"
    envmodules: *get_params("jellyfish_histo", "envmodules")
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/histo"


rule jellyfish_kmer_count_pairs:
    output:
        tsv = "{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.tsv"
    input:
        assembly = "{results}/jellyfish/{assembly}/merged.{kmer}mer_counts.jf",
        reads = "{results}/jellyfish/{analysis}/merged.{kmer}mer_counts.jf"
    conda:
        "../envs/jellyfish-kmer-utils.yaml"
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_kmer_count_pairs", "runtime", attempt)
    threads: get_params("jellyfish_kmer_count_pairs", "threads")
    log: "logs/{results}/jellyfish/{assembly}.{analysis}.{kmer}_jf.log"
    shell:
        "kmer_count_pairs {input.assembly} {input.reads} > {output.tsv}"


rule jellyfish_kmer_pairs_plot:
    """Plot kmer assembly and read pairs"""
    output:
        png = report("{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.png",
                     caption="../report/kmer_comparison.rst", category="Kmer pairs plot"),
        pngfill = report("{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.fill.png",
                         caption="../report/kmer_comparison.rst", category="Kmer pairs plot")
    input:
        tsv = "{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.tsv"
    conda:
        "../envs/jellyfish-R.yaml"
    envmodules: *get_params("jellyfish_kmer_pairs_plot", "envmodules")
    threads: get_params("jellyfish_kmer_pairs_plot", "threads")
    log: "logs/{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.png.log"
    script:
        "../scripts/assemblyeval_jellyfish_kmer_pairs_plot.R"
