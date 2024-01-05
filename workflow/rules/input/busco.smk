def all_busco_input(wildcards):
    val = []
    for analysis in cfg.analyses_w_tool("busco"):
        val.extend(analysis.tools["busco"].input)
    return val


def busco_lineage_input(wildcards, input):
    data = pd.read_table(str(input), header=None, index_col=0)
    date = data.loc[wildcards.lineage, 1]
    url = f"https://busco-data.ezlab.org/v5/data/lineages/{wildcards.lineage}.{date}.tar.gz"
    return url
