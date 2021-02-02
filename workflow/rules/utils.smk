##############################
## Map data sources
##############################
datasource_files = []
source_files = []
if datasources is not None:
    datasource_files = datasources["data"].tolist()
    source_files = datasources["source"].tolist()
src_map = dict(zip(datasource_files, source_files))
cmd_map = {'': 'ln -s', 'rsync': 'rsync -av', 'sftp': 'cp', 'file': 'ln -s'}


rule assemblyeval_get_external:
    """Get external resources with appropriate method"""
    wildcard_constraints:
        resource_file = "({})".format("|".join(x for x in datasource_files))
    output:
        "{resource_file}"
    input: uri = lambda wildcards: assemblyeval_get_external_input(src_map[wildcards.resource_file])
    params:
        cmd = lambda wildcards: cmd_map[get_uri_scheme(src_map[wildcards.resource_file])],
        uri = lambda wildcards: parse_uri(src_map[wildcards.resource_file])
    log: "logs/{resource_file}.log"
    shell:
        "{params.cmd} {params.uri} {output}"


rule assemblyeval_samtools_faidx:
    """Run samtools faidx on fasta file"""
    output:
        "{prefix}{fa}{gz}.fai"
    input:
        "{prefix}{fa}{gz}"
    resources:
        runtime = lambda wildcards, attempt: resources("assemblyeval_samtools_faidx", "runtime", attempt)
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{prefix}{fa}{gz}.fai.log"
    threads:
        1
    wrapper:
        f"{SMK_WRAPPER_PREFIX}/bio/samtools/faidx"



rule assemblyeval_save_config:
    """Save assemblyeval configuration"""
    output: report("config/assemblyeval.config.yaml", caption="../report/config.rst", category="Configuration")
    log: "logs/assemblyeval/assemblyeval_save_config.log"
    script: "../scripts/assemblyeval_save_config.py"


localrules: assemblyeval_get_external
