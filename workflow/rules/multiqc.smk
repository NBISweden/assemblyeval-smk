rule multiqc:
    output: "{reports}/multiqc.html"
    input: unpack(all_multiqc)
    resources:
        runtime = lambda wildcards, attempt: resources("multiqc", "runtime", attempt)
    params: get_params("multiqc", "options")
    log: "logs/{reports}/multiqc.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/multiqc"
