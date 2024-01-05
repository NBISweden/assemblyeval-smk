rule gmap_build:
    """Create gmap database from sequence file"""
    output:
        db=__INTERIM__ / "gmap/db/{assembly}.db.ok",
    input:
        seq=lambda wildcards: cfg.get_assembly(wildcards.assembly),
    cache: True
    threads: cfg.ruleconf("gmap_build").xthreads
    log:
        "logs/gmap_build/{assembly}.log",
    envmodules:
        *cfg.ruleconf("gmap_build").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/gmap/build")


rule gmap_map:
    """Map transcriptome to sequence database"""
    output:
        res=report(
            __INTERIM__ / "gmap/map/{assembly}-{transcriptome}.psl",
            caption="../report/gmap.rst",
            category="Gmap mapping",
        ),
    input:
        db=__INTERIM__ / "gmap/db/{assembly}.db.ok",
        transcriptome=lambda wildcards: cfg.get_transcriptome(wildcards.transcriptome),
        log="logs/gmap_build/{assembly}.log",
    threads: cfg.ruleconf("gmap_map").xthreads
    log:
        "logs/gmap_map/{assembly}-{transcriptome}.psl.log",
    envmodules:
        *cfg.ruleconf("gmap_map").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/gmap/map")
