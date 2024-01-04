def all_busco_input(wildcards):
    val = []
    for analysis in cfg.analyses_w_tool("busco"):
        val.extend(analysis.tools["busco"].input)
    return val
