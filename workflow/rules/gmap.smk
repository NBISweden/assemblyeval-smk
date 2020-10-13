rule gmap_build:
    """Create gmap database from sequence file"""
    output:
        db = __INTERIM__ / "gmap/db/{assembly}.db.ok"
    input:
        seq = get_assembly
    conda:
        "../envs/gmap.yaml"
    cache: True
    resources:
        runtime = lambda wildcards, attempt: attempt * config["gmap"]["build"]["runtime"],
        mem_mb = config["gmap"]["build"]["mem_mb"]
    threads:
        1
    log:
        "logs/gmap_build/{assembly}.log"
    wrapper:
        f"{WRAPPER_PREFIX}/bio/gmap/build"

rule gmap_map:
    """Map transcriptome to sequence database"""
    output:
        res = report(__INTERIM__ / "gmap/map/{assembly}-{transcriptome}.psl", caption="../report/gmap.rst", category="Gmap mapping")
    input:
        db = __INTERIM__ / "gmap/db/{assembly}.db.ok",
        transcriptome = get_transcriptome,
        log = "logs/gmap_build/{assembly}.log"
    conda:
        "../envs/gmap.yaml"
    resources:
        runtime = lambda wildcards, attempt: attempt * config["gmap"]["map"]["runtime"],
        mem_mb = config["gmap"]["map"]["mem_mb"]
    threads:
        lambda wildcards, attempt: attempt * config["gmap"]["map"]["threads"]
    log:
        "logs/gmap_map/{assembly}-{transcriptome}.psl.log"
    wrapper:
        f"{WRAPPER_PREFIX}/bio/gmap/map"
