rule jellyfish_count:
    output: jf = temp("{results}/qc/jellyfish/{assembly}{gz}.{kmer}mer_counts.jf")
    input: seq = get_assembly
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_count", "runtime", attempt)
    params:
        options = get_params("jellyfish_count", "options")
    threads: get_params("jellyfish_count", "threads")
    log: "logs/{results}/qc/jellyfish/{assembly}{gz}.{kmer}mer_counts.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/count"


rule jellyfish_histo:
    output: hist = "{prefix}.{kmer}_jf.hist"
    input: counts = "{prefix}.{kmer}mer_counts.jf"
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_histo", "runtime", attempt)
    params:
        options = get_params("jellyfish_histo", "options")
    threads: get_params("jellyfish_histo", "threads")
    log: "logs/{prefix}.{kmer}mer_counts.jf.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/histo"
