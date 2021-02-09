rule kraken2_parallel:
    """Run kraken2 in parallel"""
    output:
        output = temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.output.txt.gz"),
        unclassified = temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.unclassified.fasta.gz"),
        report = temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.report.txt"),
    input:
        "{interim}/kraken2/{assembly}/{db}.{length}.{partition}.fasta"
    params:
        db = get_workflow_params("kraken2", "db"),
        options = get_params("kraken2_parallel", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("kraken2_parallel", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("kraken2_parallel", "mem_mb", attempt),
    threads:
        get_params("kraken2_parallel", "threads")
    log: "logs/{interim}/kraken2/{assembly}/{db}.{length}.{partition}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/kraken2/kraken2"


rule kraken2_bedtools_make_windows:
    """Convert assembly fai to bed file defining windows"""
    output:
        bed = temp("{interim}/kraken2/{assembly}/{db}.{length}.bed")
    input:
        seq = get_assembly_index
    log: "logs/{interim}/kraken2/{assembly}/{db}.{length}.bed.log"
    shell:
        "cat {input.seq} | awk -v OFS='\\t' '{{print $1, $2}}' | bedtools makewindows -g - -w {wildcards.length} > {output.bed}"


rule kraken2_python_make_windows:
    """Partition input sequence file into windows in batches"""
    output:
        temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.fasta")
    input:
        fasta = get_assembly,
        bed = "{interim}/kraken2/{assembly}/{db}.{length}.bed"
    params:
        npart = get_workflow_params("kraken2", "npartitions")
    conda:
        "../envs/pybedtools.yaml"
    log: "logs/{interim}/kraken2/{assembly}/{db}.{length}.{partition}.fasta.log"
    script: "../scripts/assemblyeval_kraken2_python_make_windows.py"


rule kraken2_gather_results:
    output:
        output = __RESULTS__ / "kraken2/{assembly}.{db}.{length}.output.txt.gz",
        unclassified = __RESULTS__ / "kraken2/{assembly}.{db}.{length}.unclassified.fasta.gz"
    input:
        unpack(kraken2_gather_results_input)
    resources:
        runtime = lambda wildcards, attempt: resources("kraken2_gather_results", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("kraken2_gather_results", "mem_mb", attempt)
    log: "logs/kraken2/{assembly}/{db}.{length}.results.log"
    threads:
        get_params("kraken2_gather_results", "threads")
    shell:
        "cat {input.output} > {output.output}; "
        "cat {input.unclassified} > {output.unclassified}"


rule kraken2_gather_reports:
    output:
        txt = __RESULTS__ / "kraken2/{assembly}.{db}.{length}.report.txt"
    input:
        output = __RESULTS__ / "kraken2/{assembly}.{db}.{length}.output.txt.gz",
        unclassified = __RESULTS__ / "kraken2/{assembly}.{db}.{length}.unclassified.fasta.gz",
        txt = kraken2_gather_reports_input
    resources:
        runtime = lambda wildcards, attempt: resources("kraken2_gather_reports", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("kraken2_gather_reports", "mem_mb", attempt)
    log: "logs/kraken2/{assembly}/{db}.{length}.report.log"
    threads:
        get_params("kraken2_gather_reports", "threads")
    script: "../scripts/assemblyeval_kraken2_gather_reports.py"


localrules: kraken2_bedtools_make_windows, kraken2_python_make_windows
