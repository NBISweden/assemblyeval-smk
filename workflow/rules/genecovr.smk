wildcard_constraints:
    genecovr_results = str(__RESULTS__ / "genecovr")

__GENECOVR_MINMATCH__ = [0.75, 0.85, 0.9, 0.95]
__GENECOVR_NCONTIGS__ = ["bar"]
__GENECOVR_MATCH_INDEL__ = ["violin", "boxplot", "boxplot.log10"]
__GENECOVR_FN__ = ["width_violin.pdf", "qnuminsert.pdf"]
__GENECOVR_DEPTH_BREADTH__ = ["coverage", "jitter", "hist", "seqlengths"]
__GENECOVR_CSV_GZ__ = ["gene_body_coverage.csv.gz", "psldata.csv.gz", "gbc_summary.csv.gz", "ncontigs_per_transcripts.csv.gz"]


rule all_genecovr:
    input: all_genecovr_input


rule genecovr_make_csv_inputfile:
    """Generate csv input file from dataset key"""
    wildcard_constraints:
        outprefix = "({})".format("|".join([config["genecovr"][x]["outprefix"] for x in config["genecovr"].keys() if x.startswith("dataset")]))
    output:
        csv = "{outprefix}.{dataset}.csv"
    input:
        unpack(genecovr_make_csv_inputfile_input)
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{outprefix}.{dataset}.log"
    script: "../scripts/assemblyeval_genecovr_make_csv_inputfile_input.py"


rule genecovr_run:
    """Calculate gene body coverage"""
    output:
        genecovr_output()
    input:
        unpack(get_genecovr_input)
    params:
        options = get_params("genecovr_run", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("genecovr_run", "runtime", attempt)
    log: "logs/{genecovr_results}/{dataset}.log"
    threads:
        lambda wildcards, attempt: resources("genecovr_run", "threads", attempt)
    wrapper:
        f"{WRAPPER_PREFIX}/bio/genecovr"


localrules: genecovr_make_csv_inputfile
