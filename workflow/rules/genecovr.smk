rule all_genecovr:
    input:
        all_genecovr,


rule genecovr_run:
    """Calculate gene body coverage"""
    output:
        report(
            expand(
                "{{results}}/genecovr/{{analysis}}/gene_body_coverage.minmatch.{mm}.pdf",
                mm=Genecovr.__GENECOVR_MINMATCH__,
            ),
            caption="../report/genecovr_gbc.rst",
            category="Gene body coverages",
        ),
        report(
            expand(
                "{{results}}/genecovr/{{analysis}}/ncontigs_per_transcripts.{type}.mm0.75.pdf",
                type=Genecovr.__GENECOVR_NCONTIGS__,
            ),
            caption="../report/genecovr_ncontigs.rst",
            category="Number of contigs per transcript",
        ),
        report(
            expand(
                "{{results}}/genecovr/{{analysis}}/depth_breadth_{type}.mm0.75.pdf",
                type=Genecovr.__GENECOVR_DEPTH_BREADTH__,
            ),
            caption="../report/genecovr_depth_breadth.rst",
            category="Depth and breadth of coverage",
        ),
        report(
            expand(
                "{{results}}/genecovr/{{analysis}}/match_indel.{type}.pdf",
                type=Genecovr.__GENECOVR_MATCH_INDEL__,
            ),
            caption="../report/genecovr_match_indel.rst",
            category="Match and indel distributions",
        ),
        report(
            expand(
                "{{results}}/genecovr/{{analysis}}/{fn}", fn=Genecovr.__GENECOVR_FN__
            ),
            caption="../report/genecovr_match_indel.rst",
            category="Match and indel distributions",
        ),
        report(
            expand(
                "{{results}}/genecovr/{{analysis}}/{fn}",
                fn=Genecovr.__GENECOVR_CSV_GZ__,
            ),
            caption="../report/genecovr_data.rst",
            category="Data files",
        ),
    input:
        unpack(get_genecovr_input),
    params:
        options=cfg.ruleconf("genecovr_run").options,
    resources:
        runtime=cfg.ruleconf("genecovr_run").xruntime,
    log:
        "logs/{results}/genecovr/{analysis}.log",
    threads: cfg.ruleconf("genecovr_run").xthreads
    wrapper:
        f"{WRAPPER_PREFIX}/bio/genecovr"
