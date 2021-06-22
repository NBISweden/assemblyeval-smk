def all_jellyfish(wildcards):
    val = []
    for analysis in cfg.analyses_w_tool("jellyfish"):
        jellyfish = analysis.tools["jellyfish"]
        val.extend(jellyfish.input)
    return val


def jellyfish_count_input(wildcards):
    return cfg.analysis(f"{wildcards.analysis}").tools["jellyfish"].count_input(wildcards)


def jellyfish_merge_input(wildcards):
    return cfg.analysis(f"{wildcards.analysis}").tools["jellyfish"].merge_input(wildcards)
