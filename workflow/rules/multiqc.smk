rule multiqc:
    output:
        "{reports}/multiqc.html",
    input:
        unpack(all_multiqc),
    params:
        cfg.ruleconf("multiqc").options,
    log:
        "logs/{reports}/multiqc.log",
    envmodules:
        *cfg.ruleconf("multiqc").envmodules,
    wrapper:
        os.path.join(WRAPPER_PREFIX, "bio/multiqc")
