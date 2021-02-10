rule gmap_build:
    """Create gmap database from sequence file"""
    output:
        db = __INTERIM__ / "gmap/db/{assembly}.db.ok"
    input:
        seq = get_assembly
    cache: True
    resources:
        runtime = lambda wildcards, attempt: resources("gmap_build", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("gmap_build", "mem_mb", attempt),
    threads:
        lambda wildcards, attempt: resources("gmap_build", "threads", attempt)
    log:
        "logs/gmap_build/{assembly}.log"
    envmodules: *get_params("gmap_build", "envmodules")
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
    resources:
        runtime = lambda wildcards, attempt: resources("gmap_map", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("gmap_map", "mem_mb", attempt),
    threads:
        lambda wildcards, attempt: resources("gmap_map", "threads", attempt)
    log:
        "logs/gmap_map/{assembly}-{transcriptome}.psl.log"
    envmodules: *get_params("gmap_map", "envmodules")
    wrapper:
        f"{WRAPPER_PREFIX}/bio/gmap/map"
