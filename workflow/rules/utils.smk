rule assemblyeval_samtools_faidx:
    """Run samtools faidx on fasta file"""
    output:
        "{prefix}{fa}{gz}.fai",
    input:
        "{prefix}{fa}{gz}",
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{prefix}{fa}{gz}.fai.log",
    threads: 1
    wrapper:
        os.path.join(SMK_WRAPPER_PREFIX, "bio/samtools/faidx")


rule assemblyeval_save_config:
    """Save assemblyeval configuration"""
    output:
        report(
            "config/assemblyeval.config.yaml",
            caption="../report/config.rst",
            category="Configuration",
        ),
    log:
        "logs/assemblyeval/assemblyeval_save_config.log",
    script:
        "../scripts/assemblyeval_save_config.py"
