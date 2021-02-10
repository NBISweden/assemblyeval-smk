rule multiqc:
    output: "{reports}/multiqc.html"
    input: unpack(all_multiqc)
    resources:
        runtime = lambda wildcards, attempt: resources("multiqc", "runtime", attempt)
    params: get_params("multiqc", "options")
    log: "logs/{reports}/multiqc.log"
    envmodules: *get_params("multiqc", "envmodules")
    wrapper: f"{WRAPPER_PREFIX}/bio/multiqc"
