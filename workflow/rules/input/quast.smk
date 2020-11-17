def all_quast(wildcards):
    pfx = str(__RESULTS__ / "qc/quast/{assembly}/{rpt}")
    rpt = ["report.tsv", "transposed_report.tsv", "report.txt",
           "transposed_report.txt"]
    val = {
        'tsv': expand(pfx, assembly=make_assembly_ids(), rpt="report.tsv"),
        'tsv.trans': expand(pfx, assembly=make_assembly_ids(), rpt="transposed_report.tsv"),
        'txt': expand(pfx, assembly=make_assembly_ids(), rpt="report.txt"),
        'txt.trans': expand(pfx, assembly=make_assembly_ids(), rpt="report.txt")
    }
    return val
