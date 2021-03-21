def all_quast(wildcards):
    tsv = tsv_trans = txt = txt_trans = []
    for analysis, toolconf, assembly_ids, _, _ in iter_analyses("quast"):
        pfx = __RESULTS__ / f"quast/{analysis}/{{assembly}}/{{rpt}}"
        tsv.extend(expand(pfx, assembly=assembly_ids, rpt="report.tsv"))
        tsv_trans.extend(expand(pfx, assembly=assembly_ids, rpt="transposed_report.tsv"))
        txt.extend(expand(pfx, assembly=assembly_ids, rpt="report.txt"))
        txt_trans.extend(expand(pfx, assembly=assembly_ids, rpt="report.txt"))
    val = {
        'tsv': tsv,
        'tsv.trans': tsv_trans,
        'txt': txt,
        'txt.trans': txt_trans
    }
    return val
