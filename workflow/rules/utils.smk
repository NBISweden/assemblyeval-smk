##############################
## Map resource file to source file
##############################
resource_files = assemblies[assemblies["source"].notnull()]["fasta"].tolist() + \
    transcripts[transcripts["source"].notnull()]["fasta"].tolist() + \
    reads[reads["source"].notnull()]["reads"].tolist()
source_files = assemblies[assemblies["source"].notnull()]["source"].tolist() + \
    transcripts[transcripts["source"].notnull()]["source"].tolist() + \
    reads[reads["source"].notnull()]["source"].tolist()
src_map = dict(zip(resource_files, source_files))

## FIXME: add support for uris, e.g. ensembl
rule assemblyeval_get_external:
    """Get external resources with appropriate method"""
    wildcard_constraints:
        resource_file = "({})".format("|".join(x for x in resource_files))
    output:
        "{resource_file}"
    input: lambda wildcards: src_map[wildcards.resource_file]
    shell:
        "rsync -av {input} {output}"


rule assemblyeval_samtools_faidx:
    """Run samtools faidx on fasta file"""
    output:
        "{prefix}{fa}{gz}.fai"
    input:
        "{prefix}{fa}{gz}"
    resources:
        runtime = lambda wildcards, attempt: attempt * config["samtools"]["faidx"]["runtime"]
    conda:
        "../envs/samtools.yaml"
    threads:
        1
    wrapper:
        str(SMK_WRAPPER_PREFIX / "bio/samtools/faidx")


localrules: assemblyeval_get_external
