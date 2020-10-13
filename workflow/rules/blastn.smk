def _blast_blastn_all(wildcards):
    return []
    # [str(__BTK_PATH__ / blobdir / Path(x).name.rstrip(".fasta.gz")) + config["btk"]["blast"]["db"]]

rule blast_blastn_all:
    input: _blast_blastn_all


rule blast_blastn_blastdb:
    """Run blastn to a database"""
    output:
        out = "{interim}/btk/{blobdir}/{prefix}.{db}.out.gz"
    input:
        fasta = "{interim}/btk/{blobdir}/{prefix}.fasta.gz",
        blobdir = "{interim}/btk/{blobdir}/.ok"
    wildcard_constraints:
        db = config["btk"]["blast"]["db"],
        prefix = "[^/]+"
    params:
        outfmt = config["btk"]["blast"]["outfmt"],
        options = config["btk"]["blast"]["options"]
    resources:
        runtime = config["btk"]["blast"]["runtime"]
    conda:
        "blast.yaml"
    envmodules:
        config["btk"]["blast"]["envmodules"]
    threads:
        config["btk"]["blast"]["threads"]
    shell:
        "blastn {params.options} -db {wildcards.db} -query <(gzip -d -c {input.fasta}) -outfmt \"{params.outfmt}\" -num_threads {threads} | gzip -v > {output.out}"


rule blast_blastn_blastdb_parallel:
    """Run blastn in parallel to a database"""
    output:
        out = "{interim}/btk/{blobdir}/blastn/{region}/{prefix}.fasta.bgz.{db}.out.gz"
    input:
        fasta = "{interim}/btk/{blobdir}/blastn/{region}/{prefix}.fasta.bgz"
    wildcard_constraints:
        db = config["btk"]["blast"]["db"],
        prefix = "[^/]+"
    params:
        outfmt = config["btk"]["blast"]["outfmt"],
        options = config["btk"]["blast"]["options"]
    resources:
        runtime = config["btk"]["blast"]["runtime"]
    conda:
        "blast.yaml"
    envmodules:
        config["btk"]["blast"]["envmodules"]
    threads:
        config["btk"]["blast"]["threads"]
    shell:
        "blastn {params.options} -db {wildcards.db} -query <(bgzip -d -c {input.fasta}) -outfmt \"{params.outfmt}\" -num_threads {threads} | gzip -v > {output.out}"

rule seqtk_select_regions_bgzip:
    """Subselect regions"""
    output:
        bgz = temp("{interim}/btk/{blobdir}/blastn/{region}/{prefix}.fasta.bgz"),
        gzi = temp("{interim}/btk/{blobdir}/blastn/{region}/{prefix}.fasta.bgz.gzi")
    input:
        done = "{interim}/btk/{blobdir}/blastn/{prefix}.fasta.bgz_partition_bed.done",
        fasta = "{interim}/btk/{blobdir}/{prefix}.fasta.bgz"
    params:
        bed = lambda wildcards, input: os.path.join(os.path.dirname(input.done), wildcards.region, wildcards.prefix + ".fasta.bgz.bed")
    resources:
        runtime = config["btk"]["seqtk"]["runtime"]
    wildcard_constraints:
        prefix = "[^/]+",
    conda:
        "seqtk.yaml"
    threads:
        1
    priority:
        10
    shell:
        "seqtk subseq {input.fasta} {params.bed} | sed \"s/:.*//g\" | bgzip -i -I {output.gzi} -c > {output.bgz}"


def _blastn_blast_partitions(wildcards):
    print(wildcards)
    return []

rule blastn_blast_partitions:
    """Pseudo rule to run blast on a subset of partitions"""
    output:
        out = "{interim}/btk/{blobdir}/{prefix}.fasta.bgz.{db}.merged.out.gz"
    input:
        #out = _blastn_blast_partitions
        out = expand("{{interim}}/btk/{{blobdir}}/blastn/{region}/{{prefix}}.fasta.bgz.{{db}}.out.gz", region=[1+p for p in range(config["btk"]["blast"]["npartitions"])])
    wildcard_constraints:
        prefix = "[^/]+"
    resources:
        runtime = lambda wildcards, attempt: attempt * config["btk"]["blast"]["runtime"]
    run:
        import shutil
        with open(output.out, 'wb') as wfp:
            for fn in input.out:
                with open(fn, 'rb') as rfp:
                    shutil.copyfileobj(rfp, wfp)

rule python_partition_bed_file:
    """Partition bed file greedily into a subset of files"""
    output:
        #target = temp("{interim}/{assembly}/blastn/{partition}/{prefix}.fa.bgz.bed")
        done = "{interim}/btk/{blobdir}/blastn/{prefix}.fasta.bgz_partition_bed.done"
    input:
        targets = "{interim}/btk/{blobdir}/{prefix}.fasta.bgz.bed",
        genome_size = "{interim}/btk/{blobdir}/{prefix}.fasta.bgz.genome_size.txt",
    params:
        npart = config["btk"]["blast"]["npartitions"]
    resources:
        runtime = 60
    priority:
        50
    wildcard_constraints:
        prefix = "[^/]+"
    threads:
        1
    conda:
        "pybedtools.yaml"
    script:
        "../../scripts/partition_bed.py"


rule unix_create_genome_size_file:
    """Create bedtools genome file"""
    output:
        "{interim}/btk/{blobdir}/{prefix}.fasta.gz.genome_size.txt"
    input:
        "{interim}/btk/{blobdir}/{prefix}.fasta.bgz.fai"
    wildcard_constraints:
        prefix = "[^/]+"
    conda:
        "pybedtools.yaml"
    shell:
        "cat {input} | awk '{{OFS=\"\t\"; print $1, $2}}' > {output}"
