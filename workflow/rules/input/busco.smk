def all_busco_input(wildcards):
    val = []
    for analysis, toolconf, assembly_ids, _, _ in iter_analyses("busco"):
        pfx = __RESULTS__ / f"busco/{analysis}/{{assembly}}/{{mode}}/run_{{lineage}}/short_summary_{{assembly}}.txt"
        tmp = expand(pfx, assembly=assembly_ids, mode=toolconf.get("mode", "genome"),
                     lineage=toolconf["lineage"])
        val.extend(tmp)
    return val
