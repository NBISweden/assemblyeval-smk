def get_btk_all(wildcards):
    retval = []
    for v in config["btk"].keys():
        if not v.startswith("blobdir"):
            continue
        if len(config["btk"][v]["fasta"]) > 0:
            retval.append(str(__INTERIM__ / f"btk/{v}/gc.json"))
        for bam in config["btk"][v]["bam"]:
            retval.append(str(__INTERIM__ / f"btk/{v}/{Path(bam).name}.sort_cov.json"))
        if len(config["btk"][v]["bls"]) > 0:
            retval.append(str(__INTERIM__ / f"btk/{v}/bestsumorder_class_score.json"))
    return retval[0]


def _btk_link_fasta_input(wildcards):
    return [_btk_link_fasta_input_paths(wildcards)["path"]]


def _btk_link_fasta_input_paths(wildcards):
    ret = {}
    fn = (
        Path(wildcards.interim)
        / "btk"
        / wildcards.blobdir
        / f"{wildcards.prefix}.fasta.gz"
    )
    for x in config["btk"][wildcards.blobdir]["fasta"]:
        if str(fn.name) == str(Path(x).name):
            ret["path"] = str(Path(x))
            ret["abspath"] = str(Path(x).absolute())
            break
    return ret
