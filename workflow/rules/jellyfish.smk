rule jellyfish_count:
    output: jf = "{results}/jellyfish/{datatype}/{prefix}{gz}.{kmer}mer_counts.jf"
    input: unpack(jellyfish_count_input)
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_count", "runtime", attempt)
    params:
        options = get_params("jellyfish_count", "options"),
        tmpdir = get_workflow_params("jellyfish", "tmpdir")
    wildcard_constraints:
        datatype = "(assembly|reads)"
    threads: get_params("jellyfish_count", "threads")
    log: "logs/{results}/jellyfish/{datatype}/{prefix}{gz}.{kmer}mer_counts.log"
    envmodules: *get_params("jellyfish_count", "envmodules")
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/count"


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


rule jellyfish_plot:
    output:
        png = report("{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.png",
                     caption="../report/kmer_comparison.rst", category="Kmer comparison"),
        tsv = "{results}/jellyfish/kmer_comparison/{assembly}.{analysis}.{kmer}_jf.tsv"
    input:
        assembly = "{results}/jellyfish/assembly/{assembly}.{kmer}mer_counts.jf",
        reads = "{results}/jellyfish/reads/{analysis}.{kmer}mer_counts.jf"
    params:
        tmpdir = get_workflow_params("jellyfish", "tmpdir")
    conda:
        "../envs/jellyfish-python.yaml"
    resources:
        runtime = lambda wildcards, attempt: resources("jellyfish_plot", "runtime", attempt)
    threads: get_params("jellyfish_plot", "threads")
    log: "logs/{results}/jellyfish/{assembly}.{analysis}.{kmer}_jf.log"
    script:
        "../scripts/assemblyeval_jellyfish_plot.py"
