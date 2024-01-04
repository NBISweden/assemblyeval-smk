# include: "blastn.smk"
# Rule to collect all btk for one blobdir
# rule btk_all:
#     """btk collect all outputs"""
#     input: get_btk_all
# rule btk_create_blobdir:
#     """btk create blobdir."""
#     output:
#         flag = "{interim}/btk/{blobdir}/.ok"
#     params:
#         blobtools = config["btk"]["exe"]
#     conda:
#         "../envs/btk.yaml"
#     log: "logs/{interim}/btk/{blobdir}.ok.log"
#     container:
#         "docker://genomehubs/blobtoolkit"
#     shell:
#         "{params.blobtools} create $(dirname {output.flag}) && touch {output.flag}"
# rule btk_link_fasta:
#     """btk link fasta file from config"""
#     params:
#         abspath = lambda wildcards: _btk_link_fasta_input_paths(wildcards)['abspath']
#     output:
#         "{interim}/btk/{blobdir}/{prefix}.fasta.gz"
#     input:
#         _btk_link_fasta_input
#     log: "logs/{interim}/btk/{blobdir}/{prefix}.fasta.gz.log"
#     wildcard_constraints:
#         prefix = "[^/]+"
#     shell:
#         "ln -s {params.abspath} {output}"
# rule btk_convert_fasta_gz_to_bgzip:
#     """Convert fasta.gz to fasta.bgz"""
#     output:
#         bgz = "{interim}/btk/{blobdir}/{prefix}.fasta.bgz",
#         gzi = "{interim}/btk/{blobdir}/{prefix}.fasta.bgz.gzi"
#     input:
#         gz = "{interim}/btk/{blobdir}/{prefix}.fasta.gz"
#     conda:
#         "blast.yaml"
#     log: "logs/{interim}/btk/{blobdir}/{prefix}.fasta.bgz.log"
#     wildcard_constraints:
#         prefix = "[^/]+"
#     threads:
#         1
#     resources:
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["convert_fasta_gz_to_bgzip"]["runtime"]
#     shell:
#         "zcat {input.gz} | bgzip -i -I {output.gzi} -c > {output.bgz}"
# rule btk_add_fasta:
#     """btk add fasta"""
#     output:
#         expand("{{interim}}/btk/{{blobdir}}/{out}", out=["gc.json", "identifiers.json", "length.json", "meta.json", "ncount.json"])
#     input:
#         fasta = lambda wildcards: config["btk"][wildcards.blobdir]["fasta"],
#         blobdir = "{interim}/btk/{blobdir}/.ok"
#     params:
#         blobtools = config["btk"]["exe"]
#     resources:
#         mem_mb = config["btk"]["add"]["mem_mb"],
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["add_fasta"]["runtime"]
#     log: "logs/{interim}/btk/{blobdir}/add_fasta.log"
#     conda:
#         "../envs/btk.yaml"
#     container:
#         "docker://genomehubs/blobtoolkit"
#     shell:
#         "{params.blobtools} add --replace --fasta {input.fasta} $(dirname {input.blobdir})"
# rule btk_samtools_index:
#     """Index bam file"""
#     output:
#         bai = "{interim}/btk/{blobdir}/{prefix}.bai"
#     input:
#         bam = "{interim}/btk/{blobdir}/{prefix}.bam"
#     conda:
#         "blast.yaml"
#     resources:
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["index"]["runtime"]
#     threads:
#         config["btk"]["index"]["threads"]
#     shell: "samtools index -@ {threads} -b {input.bam} {output.bai}"
# rule btk_samtools_sort:
#     """Sort input file"""
#     output:
#         bam = temp("{interim}/btk/{blobdir}/{prefix}.sort.bam")
#     input:
#         bam = "{interim}/btk/{blobdir}/{prefix}.bam"
#     params:
#         options = config["btk"]["sort"]["options"]
#     resources:
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["sort"]["runtime"]
#     threads:
#         config["btk"]["sort"]["threads"]
#     conda:
#         "blast.yaml"
#     shell: "samtools sort {params.options} -@ {threads} {input.bam} -o {output.bam}"
# rule btk_samtools_faidx:
#     """samtools fasta index fasta file"""
#     output:
#         fai = "{interim}/btk/{blobdir}/{prefix}.fasta.{gz}.fai"
#     input:
#         fa = "{interim}/btk/{blobdir}/{prefix}.fasta.{gz}"
#     conda:
#         "blast.yaml"
#     resources:
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["faidx"]["runtime"]
#     wildcard_constraints:
#         gz = "(gz|bgz)"
#     threads:
#         config["btk"]["faidx"]["threads"]
#     shell: "samtools faidx {input.fa} -o {output.fai}"
# rule btk_create_genome_size_file:
#     """Create bedtools genome size file"""
#     output:
#         "{interim}/btk/{blobdir}/{prefix}.fasta.{gz}.genome_size.txt"
#     input:
#         "{interim}/btk/{blobdir}/{prefix}.fasta.{gz}.fai"
#     wildcard_constraints:
#         gz = "(gz|bgz)",
#         prefix = "[^/]+"
#     shell:
#         "cat {input} | awk '{{OFS=\"\t\"; print $1, $2}}' > {output}"
# rule btk_unix_awk_genome_size_to_bed:
#     """Convert genome_size to bed"""
#     output:
#         bed = "{interim}/btk/{prefix}.bed"
#     input:
#         gs = "{interim}/btk/{prefix}.genome_size.txt"
#     shell:
#         "cat {input.gs} | awk '{{OFS=\"\t\"; print $1, 0, $2}}' > {output.bed}"
# rule btk_link_mapped_reads:
#     output:
#         bam = "{interim}/btk/{blobdir}/{prefix}.bam"
#     input:
#         bam = lambda wildcards: [str(Path(x)) for x in config["btk"][wildcards.blobdir]["bam"] if Path(x).name == f"{wildcards.prefix}.bam"][0]
#     params:
#         bam = lambda wildcards, input: str(Path(input.bam).absolute())
#     shell:
#         "ln -s {params.bam} {output.bam}"
# rule btk_add_cov_inputs:
#     """btk cov input files"""
#     input: [str(__BTK_PATH__ / blobdir / Path(x).name.rstrip(".bam")) + ".sort_cov.json" for blobdir in config["btk"].keys() if blobdir.startswith("blobdir") for x in config["btk"][blobdir]["bam"]]
# rule btk_add_cov:
#     """btk add coverage results from sorted bam files"""
#     output:
#         expand("{{interim}}/btk/{{blobdir}}/{{prefix}}.sort_{sfx}.json", sfx=["cov", "read_cov"])
#     input:
#         bam = "{interim}/btk/{blobdir}/{prefix}.sort.bam",
#         bai = "{interim}/btk/{blobdir}/{prefix}.sort.bai",
#         identifiers = "{interim}/btk/{blobdir}/identifiers.json"
#     conda:
#         "btk.yaml"
#     params:
#         blobtools = config["btk"]["exe"]
#     resources:
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["add_cov"]["runtime"],
#         mem_mb = config["btk"]["add_cov"]["mem_mb"]
#     log:
#         "{interim}/btk/{blobdir}/{prefix}.sort_cov.log"
#     threads:
#         config["btk"]["add_cov"]["threads"]
#     shell:
#         "{params.blobtools} add --replace --cov {input.bam} --threads {threads} $(dirname {output[0]}) 2>&1 > {log}"
# rule btk_add_blastn_all:
#     """btk blastn input files"""
#     input: [str(__BTK_PATH__ / blobdir / "bestsumorder_class_score.json") for blobdir in config["btk"].keys() if blobdir.startswith("blobdir") for x in config["btk"][blobdir]["bls"]]
# rule btk_add_blastn:
#     """btk add blastn results"""
#     output:
#         "{interim}/btk/{blobdir}/bestsumorder_class_score.json"
#     input:
#         identifiers = "{interim}/btk/{blobdir}/identifiers.json",
#         taxdump = config["btk"]["taxdump"],
#         bls = lambda wildcards: config["btk"][wildcards.blobdir]["bls"]
#     resources:
#         runtime = lambda wildcards, attempt: attempt * config["btk"]["add_blastn"]["runtime"]
#     threads:
#         config["btk"]["add_blastn"]["threads"]
#     params:
#         blobtools = config["btk"]["exe"],
#         hits_cols = config["btk"]["hits-cols"]
#     container:
#         "docker://genomehub/blobtoolskit"
#     conda:
#         "../envs/btk.yaml"
#     shell:
#          "{params.blobtools} add --replace --hits {input.bls} --taxrule bestsumorder --hits-cols {params.hits_cols} --taxdump $(dirname {input.taxdump}) $(dirname {output})"
# localrules: btk_link_fasta, btk_create_genome_size_file, btk_unix_awk_genome_size_to_bed, btk_link_mapped_reads, btk_create_blobdir
# ruleorder: btk_samtools_sort > btk_link_mapped_reads
