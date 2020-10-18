##############################
## Map resource file to source file
##############################
resource_files = assemblies[assemblies["source"].notnull()]["fasta"].tolist()
source_files = assemblies[assemblies["source"].notnull()]["source"].tolist()
if transcripts is not None:
    resource_files += transcripts[transcripts["source"].notnull()]["fasta"].tolist()
    source_files += transcripts[transcripts["source"].notnull()]["source"].tolist()
if reads is not None:
    resource_files += reads[reads["source"].notnull()]["reads"].tolist()
    source_files += reads[reads["source"].notnull()]["source"].tolist()
src_map = dict(zip(resource_files, source_files))
cmd_map = {'': 'ln -s', 'rsync': 'rsync -av', 'sftp': 'cp', 'file': 'ln -s'}

rule assemblyeval_get_external:
    """Get external resources with appropriate method"""
    wildcard_constraints:
        resource_file = "({})".format("|".join(x for x in resource_files))
    output:
        "{resource_file}"
    input: uri = lambda wildcards: parse_uri(src_map[wildcards.resource_file])
    params:
        cmd = lambda wildcards: cmd_map[get_uri_scheme(src_map[wildcards.resource_file])]
    log: "logs/{resource_file}.log"
    shell:
        "{params.cmd} {input} {output}"


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
    log:
        "logs/{prefix}{fa}{gz}.fai.log"
    threads:
        1
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/0.66.0/bio/samtools/faidx"


localrules: assemblyeval_get_external
