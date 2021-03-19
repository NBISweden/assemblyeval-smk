rule all_genecovr:
    input: all_genecovr


rule genecovr_run:
    """Calculate gene body coverage"""
    output:
        report(expand("{{genecovr_results}}/{{analysis}}/gene_body_coverage.minmatch.{mm}.pdf", mm=__GENECOVR_MINMATCH__),
               caption="../report/genecovr_gbc.rst", category="Gene body coverages"),
        report(expand("{{genecovr_results}}/{{analysis}}/ncontigs_per_transcripts.{type}.mm0.75.pdf",
                      type=__GENECOVR_NCONTIGS__),
               caption="../report/genecovr_ncontigs.rst", category="Number of contigs per transcript"),
        report(expand("{{genecovr_results}}/{{analysis}}/depth_breadth_{type}.mm0.75.pdf",
                      type=__GENECOVR_DEPTH_BREADTH__),
               caption="../report/genecovr_depth_breadth.rst", category="Depth and breadth of coverage"),
        report(expand("{{genecovr_results}}/{{analysis}}/match_indel.{type}.pdf",
                      type=__GENECOVR_MATCH_INDEL__),
               caption="../report/genecovr_match_indel.rst", category="Match and indel distributions"),
        report(expand("{{genecovr_results}}/{{analysis}}/{fn}", fn=__GENECOVR_FN__),
               caption="../report/genecovr_match_indel.rst", category="Match and indel distributions"),
        report(expand("{{genecovr_results}}/{{analysis}}/{fn}", fn=__GENECOVR_CSV_GZ__),
               caption="../report/genecovr_data.rst", category="Data files")
    input:
        unpack(get_genecovr_input)
    params:
        options = get_params("genecovr_run", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("genecovr_run", "runtime", attempt)
    wildcard_constraints:
        genecovr_results = str(__RESULTS__ / "genecovr")
    log: "logs/{genecovr_results}/{analysis}.log"
    threads:
        lambda wildcards, attempt: resources("genecovr_run", "threads", attempt)
    wrapper:
        f"{WRAPPER_PREFIX}/bio/genecovr"
