def all_kraken2_input(wildcards):
    val = []
    for analysis in cfg.analyses_w_tool("kraken2"):
        kraken2 = analysis.tools["kraken2"]
        val.extend(kraken2.input)
    return val


def kraken2_gather_results_input(wildcards):
    return cfg.analysis(wildcards.analysis).tools["kraken2"].results(wildcards)

def kraken2_gather_reports_input(wildcards):
    return cfg.analysis(wildcards.analysis).tools["kraken2"].reports(wildcards)
