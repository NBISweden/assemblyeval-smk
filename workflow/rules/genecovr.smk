wildcard_constraints:
    genecovr_results = str(__RESULTS__ / "genecovr")

__GENECOVR_MINMATCH__ = [0.75, 0.85, 0.9, 0.95]
__GENECOVR_NCONTIGS__ = ["bar"]
__GENECOVR_INDEL__ = ["violin", "boxplot", "boxplot.log10"]
__GENECOVR_DEPTH_BREADTH__ = ["coverage", "jitter", "hist"]
__GENECOVR_FN__ = ["gene_body_coverage.csv.gz", "psldata.csv.gz", "width_violin.pdf", "qnuminsert.pdf"]


rule genecovr_all:
    input: get_genecovr_all


def get_genecovr_input(wildcards):
    retval = []
    csvfile = config["genecovr"][f"csv_{wildcards.dataset}"]
    data = pd.read_csv(csvfile, header=None).set_index(0, drop=False)
    data.columns = ["dataset", "psl", "assembly", "trxset"]
    for k in ["psl", "assembly", "trxset"]:
        retval += [x for x in data[k].tolist() if x is not None]
    return {'csv': csvfile, 'files': retval}


rule genecovr:
    """Calculate gene body coverage"""
    output:
        expand("{{genecovr_results}}/{{dataset}}/gene_body_coverage.minmatch.{mm}.pdf", mm=__GENECOVR_MINMATCH__) + \
        expand("{{genecovr_results}}/{{dataset}}/ncontigs_per_transcripts.{type}.mm0.75.pdf", type=__GENECOVR_NCONTIGS__) + \
        expand("{{genecovr_results}}/{{dataset}}/depth_breadth_{type}.mm0.75.pdf", type=__GENECOVR_DEPTH_BREADTH__) + \
        expand("{{genecovr_results}}/{{dataset}}/match_indel.{type}.pdf", type=__GENECOVR_INDEL__) + \
        expand("{{genecovr_results}}/{{dataset}}/{fn}", fn=__GENECOVR_FN__)
    input:
        unpack(get_genecovr_input)
    params:
        options = config["genecovr"]["options"],
        exe = config["genecovr"]["exe"]
    resources:
        runtime = lambda wildcards, attempt: 2 * attempt * config["genecovr"]["runtime"]
    conda:
        "../envs/genecovr.yaml"
    threads:
        config["genecovr"]["threads"]
    shell:
        "{params.exe} -p {threads} {params.options} -d results/genecovr/{wildcards.dataset} {input.csv}"
