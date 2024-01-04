rule jellyfish_make_chunked_input:
    """Chunkify assembly input"""
    output:
        fasta=temp("{interim}/jellyfish/{analysis}/{assembly}/{partition}.fasta"),
    input:
        seq=lambda wildcards: cfg.get_assembly(wildcards.assembly),
        faidx=lambda wildcards: cfg.get_assembly(wildcards.assembly, fai=True),
    params:
        npartitions=lambda wildcards: cfg.analysis(wildcards.analysis)
        .tools["jellyfish"]
        .npartitions,
    conda:
        "../envs/pybedtools.yaml"
    resources:
        runtime=cfg.ruleconf("jellyfish_make_chunked_input").xruntime,
    threads: cfg.ruleconf("jellyfish_make_chunked_input").threads
    log:
        "logs/{interim}/jellyfish/{analysis}/{assembly}/{partition}.fasta.log",
    script:
        "../scripts/assemblyeval_pybedtools_make_chunks.py"


rule jellyfish_count_chunk:
    """Count kmers in chunk"""
    output:
        jf=temp("{interim}/jellyfish/{analysis}/{dataset}/{prefix}.{kmer}mer_counts.jf"),
    input:
        unpack(jellyfish_count_input),
    resources:
        runtime=cfg.ruleconf("jellyfish_count_chunk").xruntime,
    params:
        options=cfg.ruleconf("jellyfish_count_chunk").options,
        tmpdir=lambda wildcards: cfg.analysis(wildcards.analysis)
        .tools["jellyfish"]
        .tmpdir,
    threads: cfg.ruleconf("jellyfish_count_chunk").xthreads
    log:
        "logs/{interim}/jellyfish/{analysis}/{dataset}/{prefix}.{kmer}mer_counts.log",
    envmodules:
        *cfg.ruleconf("jellyfish_count_chunk").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/jellyfish/count")


rule jellyfish_merge:
    """Merge chunked kmer counts. Dataset refers either to assembly id or
    read id."""
    output:
        jf="{results}/jellyfish/{analysis}/{dataset}/merged.{kmer}mer_counts.jf",
    input:
        unpack(jellyfish_merge_input),
    resources:
        runtime=cfg.ruleconf("jellyfish_merge").xruntime,
    params:
        options=cfg.ruleconf("jellyfish_merge").options,
        tmpdir=lambda wildcards: cfg.analysis(wildcards.analysis)
        .tools["jellyfish"]
        .tmpdir,
        npartitions=lambda wildcards: cfg.analysis(wildcards.analysis)
        .tools["jellyfish"]
        .npartitions,
    log:
        "logs/{results}/jellyfish/{analysis}/{dataset}/merged.{kmer}mer_counts.jf.log",
    envmodules:
        *cfg.ruleconf("jellyfish_merge").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/jellyfish/merge")


rule jellyfish_histo:
    output:
        hist="{results}/jellyfish/{analysis}/{dataset}/{kmer}_jf.hist",
    input:
        counts="{results}/jellyfish/{analysis}/{dataset}/{kmer}mer_counts.jf",
    resources:
        runtime=cfg.ruleconf("jellyfish_histo").xruntime,
    params:
        options=cfg.ruleconf("jellyfish_histo").options,
        tmpdir=lambda wildcards: cfg.analysis(wildcards.analysis)
        .tools["jellyfish"]
        .tmpdir,
    threads: cfg.ruleconf("jellyfish_histo").threads
    log:
        "logs/{results}/jellyfish/{analysis}/{dataset}/{kmer}_jf.hist.log",
    envmodules:
        *cfg.ruleconf("jellyfish_histo").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/jellyfish/histo")


rule jellyfish_kmer_count_pairs:
    output:
        tsv="{results}/jellyfish/{analysis}/kmer_comparison/{assembly}.{kmer}_jf.tsv",
    input:
        assembly="{results}/jellyfish/{analysis}/{assembly}/merged.{kmer}mer_counts.jf",
        reads="{results}/jellyfish/{analysis}/kmer_comparison/merged.{kmer}mer_counts.jf",
    params:
        prefix="{results}/jellyfish/{analysis}/kmer_comparison/{assembly}.{kmer}_jf",
    conda:
        "../envs/jellyfish-kmer-utils.yaml"
    resources:
        runtime=cfg.ruleconf("jellyfish_kmer_count_pairs").xruntime,
    threads: cfg.ruleconf("jellyfish_kmer_count_pairs").threads
    log:
        "logs/{results}/jellyfish/{analysis}/{assembly}.{kmer}_jf.log",
    shell:
        "kmer_count_pairs {input.assembly} {input.reads} {params.prefix}"


rule jellyfish_kmer_pairs_plot:
    """Plot kmer assembly and read pairs"""
    output:
        png=report(
            "{results}/jellyfish/{analysis}/kmer_comparison/{assembly}.{kmer}_jf.png",
            caption="../report/kmer_comparison.rst",
            category="Kmer pairs plot",
        ),
        pngfill=report(
            "{results}/jellyfish/{analysis}/kmer_comparison/{assembly}.{kmer}_jf.fill.png",
            caption="../report/kmer_comparison.rst",
            category="Kmer pairs plot",
        ),
    input:
        tsv="{results}/jellyfish/{analysis}/kmer_comparison/{assembly}.{kmer}_jf.tsv",
    conda:
        "../envs/jellyfish-R.yaml"
    envmodules:
        *cfg.ruleconf("jellyfish_kmer_pairs_plot").envmodules,
    threads: cfg.ruleconf("jellyfish_kmer_pairs_plot").threads
    log:
        "logs/{results}/jellyfish/{analysis}/kmer_comparison/{assembly}.{kmer}_jf.png.log",
    script:
        "../scripts/assemblyeval_jellyfish_kmer_pairs_plot.R"
