rule kraken2_parallel:
    """Run kraken2 in parallel"""
    output:
        output = temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.output.txt"),
        unclassified = temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.unclassified.fasta"),
        report = temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.report.txt"),
    input:
        "{interim}/kraken2/{assembly}/{db}.{length}.{partition}.fasta"
    params:
        db = lambda wildcards: config["kraken2"]["db"],
        options = get_params("kraken2_parallel", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("kraken2_parallel", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("kraken2_parallel", "mem_mb", attempt),
    threads:
        get_params("kraken2_parallel", "threads")
    log: "logs/{interim}/kraken2/{assembly}/{db}.{length}.{partition}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/kraken2"


rule kraken2_bedtools_make_windows:
    """Convert assembly fai to bed file defining windows"""
    output:
        temp("{interim}/kraken2/{assembly}/{db}.{length}.bed")
    input:
        get_assembly_index
    log: "logs/{interim}/kraken2/{assembly}/{db}.{length}.bed.log"
    shell:
        "cat {input[0]} | awk -v OFS='\\t' '{{print $1, $2}}' | bedtools makewindows -g - -w {wildcards.length} > {output[0]}"


rule kraken2_python_make_windows:
    """Partition input sequence file into windows in batches"""
    output:
        temp("{interim}/kraken2/{assembly}/{db}.{length}.{partition}.fasta")
    input:
        fasta = get_assembly,
        bed = "{interim}/kraken2/{assembly}/{db}.{length}.bed"
    params:
        npart = config["kraken2"]["npartitions"]
    conda:
        "../envs/pybedtools.yaml"
    script: "../scripts/assemblyeval_kraken2_python_make_windows.py"


rule kraken2_gather_results:
    output:
        output = __RESULTS__ / "kraken2/{assembly}/{db}.{length}.output.txt.gz",
        unclassified = __RESULTS__ / "kraken2/{assembly}/{db}.{length}.unclassified.fasta.gz"
    input:
        output = expand(__INTERIM__ / "kraken2/{{assembly}}/{{db}}.{{length}}.{partition}.output.txt", partition=range(0, config["kraken2"]["npartitions"])),
        unclassified = expand(__INTERIM__ / "kraken2/{{assembly}}/{{db}}.{{length}}.{partition}.unclassified.fasta", partition=range(0, config["kraken2"]["npartitions"]))
    shell:
        "cat {input.output} | gzip -v - > {output.output}; "
        "cat {input.unclassified} | gzip -v - > {output.unclassified}"


rule kraken2_gather_reports:
    output:
        txt = __RESULTS__ / "kraken2/{assembly}.{db}.{length}.report.txt"
    input:
        output = __RESULTS__ / "kraken2/{assembly}/{db}.{length}.output.txt.gz",
        unclassified = __RESULTS__ / "kraken2/{assembly}/{db}.{length}.unclassified.fasta.gz",
        txt = expand(__INTERIM__ / "kraken2/{{assembly}}/{{db}}.{{length}}.{partition}.report.txt", partition=range(0, config["kraken2"]["npartitions"]))
    resources:
        runtime = lambda wildcards, attempt: resources("kraken2_gather_reports", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("kraken2_gather_reports", "mem_mb", attempt)
    log: "logs/kraken2/{assembly}/{db}.{length}.report.log"
    threads:
        get_params("kraken2_gather_reports", "threads")
    script: "../scripts/assemblyeval_kraken2_gather_reports.py"


localrules: kraken2_bedtools_make_windows, kraken2_python_make_windows, kraken2_gather_results
