wildcard_constraints:
    genecovr_results = str(__RESULTS__ / "genecovr")

__GENECOVR_MINMATCH__ = [0.75, 0.85, 0.9, 0.95]
__GENECOVR_NCONTIGS__ = ["bar"]
__GENECOVR_MATCH_INDEL__ = ["violin", "boxplot", "boxplot.log10"]
__GENECOVR_FN__ = ["width_violin.pdf", "qnuminsert.pdf"]
__GENECOVR_DEPTH_BREADTH__ = ["coverage", "dataset_jitter", "dataset_hist"]
__GENECOVR_CSV_GZ__ = ["gene_body_coverage.csv.gz", "psldata.csv.gz", "gbc_summary.csv.gz", "ncontigs_per_transcripts.csv.gz"]



rule genecovr_all:
    input: get_genecovr_all


rule genecovr_make_csv_inputfile:
    """Generate csv input file from dataset key"""
    wildcard_constraints:
        outprefix = "({})".format("|".join([config["genecovr"][x]["outprefix"] for x in config["genecovr"].keys() if x.startswith("dataset")]))
    output:
        csv = "{outprefix}.{dataset}.csv"
    input:
        unpack(genecovr_make_csv_inputfile_input)
    log:
        "logs/{outprefix}.{dataset}.log"
    run:
        df = genecovr_make_csv_inputfile_dataframe(dict(input), wildcards)
        df.to_csv(output.csv, index=False, header=False)


rule genecovr:
    """Calculate gene body coverage"""
    output:
        genecovr_output()
    input:
        unpack(get_genecovr_input)
    params:
        options = config["genecovr"]["options"],
        exe = config["genecovr"]["exe"]
    resources:
        runtime = lambda wildcards, attempt: attempt * config["genecovr"]["runtime"]
    log: "logs/{genecovr_results}/{dataset}.log"
    conda:
        "../envs/genecovr.yaml"
    threads:
        config["genecovr"]["threads"]
    shell:
        "{params.exe} -p {threads} {params.options} -d results/genecovr/{wildcards.dataset} {input.csv}"


localrules: genecovr_make_csv_inputfile
