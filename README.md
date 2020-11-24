# Snakemake workflow: Assembly evaluation workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.25.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build status](https://github.com/NBISweden/assemblyeval-smk/workflows/Tests/badge.svg?branch=main)](https://github.com/NBISweden/assemblyeval-smk/actions?query=workflow%3ATests)  ![License](https://img.shields.io/badge/license-MIT-blue.svg)

Snakemake workflow for genome assembly evaluation.

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above; currently N/A).

## Authors

* Per Unneberg (@percyfal)

## Features

- align transcripts to a reference sequence
- estimate gene body coverage with [genecovr](https://github.com/NBISweden/genecovr)
- run [quast](http://bioinf.spbau.ru/quast), [jellyfish](http://www.genome.umd.edu/jellyfish.html), and [busco](https://busco.ezlab.org)
- summarize quality metrics with [MultiQC](https://multiqc.info)
- (WIP): run some of the steps for adding reads, coverage files, and
  sequences to the [blobtoolkit](https://blobtoolkit.genomehubs.org) viewer

## Quickstart

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files
in the `config/` folder. Adjust `config.yaml` to configure the
workflow execution. In addition, edit the following files:

`assemblies.tsv`
    Assembly definition file which lists species, version, and assembly
	fasta file. Mandatory.

`transcripts.tsv`
    Transcripts definition file that lists transcript fasta files and
    connects them to a given species

`reads.tsv`
    Raw sequence read files.

`datasources.tsv`
	data-source key-value pairs definining workflow files and sources
	specified as [Uniform Resource Identifiers
	(uri)](https://en.wikipedia.org/wiki/Uniform_Resource_Identifier)

NOTE: the config directory doesn't have to be in the workflow source
directory, in which case snakemake must be invoked with the full path
to the Snakemake file:

	snakemake -s /path/to/manticore-smk/workflow/Snakefile


### Step 3: Install Snakemake

Install Snakemake using
[conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

You can also use a [snakemake
profile](https://github.com/snakemake-profiles/) for fine-tuning
executions. For instance, to use the [slurm
profile](https://github.com/Snakemake-Profiles/slurm) run

	cookiecutter https://github.com/Snakemake-Profiles/slurm.git
	snakemake --use-conda --profile slurm --jobs 100

If you not only want to fix the software stack but also the underlying OS, use

    snakemake --use-conda --use-singularity

in combination with any of the modes above. See the [Snakemake
documentation](https://snakemake.readthedocs.io/en/stable/executable.html)
for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained
interactive HTML report with all results via:

    snakemake --report report.html

The report contains documentation and results from the workflow.

### Step 6: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.

## Configuration

For a quick overview of example configuration files, see [config/config.yaml](https://github.com/NBISweden/assemblyeval-smk/blob/main/config/config.yaml) and the test configuration [.test/config/config.yaml](https://github.com/NBISweden/assemblyeval-smk/blob/main/.test/config/config.yaml)


### Schemas

All configuration files are evaluated against [configuration
schemas](https://github.com/NBISweden/assemblyeval-smk/tree/main/workflow/schemas).
The validation ensures configuration keys are populated, which
minimizes configuration overhead for the user. The schemas are
self-documented and define and describe all available configuration
options.

As an example, `workflow/schemas/assemblies.schema.yaml` defines a
tabular sample input file format for individual assemblies. There are
four defined properties: `id` (identifier), `species`, `version`, and
`fasta`, of which `species`, `version` and `fasta` are required (`id`
will be auto-generated from `species`_`version` if missing).


See the tutorial [understanding
jsonschema](https://json-schema.org/understanding-json-schema/) for an
accessible introduction to schemas.

### Workflow configuration

The workflow is configured by defining up to four input data files and
what applications to run. The input files are `assemblies`,
`transcripts`, `reads`, and `datasources`, where only `assemblies` is
mandatory:

	assemblies: config/assemblies.tsv
	transcripts: config/transcripts.tsv
	reads: config/reads.tsv
	datasources: config/datasources.tsv

Each input data file is a tab-separated file listing file names and
associated metadata. See the schemas for details or tests for
examples.

Tool configurations are listed below.

#### genecovr

The purpose of the `genecovr` property is to define what assemblies
and transcripts to use for gene body coverage calculations. Briefly,
it consists of sections that either list an input file to `genecovr`
or defines a combination of assemblies and transcripts:

	genecovr:
	  dataset1:
        csvfile: config/genecovr.csv
	  dataset2:
        assemblies: ["foo_v2", "foo_v1"]
        transcripts: ["A", "B"]

The `genecovr.csv` file consists of columns `dataset` (an identifier
name), `psl` that gives a path to a psl file of mapped transcripts,
`assembly` which points to an assembly file, and `trxset` which points
to a transcript fasta file. The `assemblies` and `transcripts` list
assembly and transcript identifiers as defined in the input assembly
and transcript files described in the previous section.

#### busco

Busco needs a lineage to run which at current is the only property
needed here.

	busco:
      lineage: viridiplantae_odb10

#### jellyfish

Jellyfish will count kmers of given sizes. The configuration section
lists what kmer sizes to use:

	jellyfish:
	  kmer: [21]

#### tools

Some tools, including `quast`, don't have any configurable options at
present. To run them, simply include the name in the `tools` property:

	tools:
	  - quast

#### btk

WIP.


### Resource configuration

Every rule has a corresponding configuration entry. Hence, it is
possible to fine tune resource configuration of `threads`, `runtime`,
and `mem_mb`, as well as modify and amend program `options`. In most
cases, these values are not initialized, in which case resources fall
back on default values defined in the `resources.default`
configuration section:

    resources.default:
      threads: 1
      mem_mb: 8192
      runtime: 120
      options: ""
      java_options: ""
      java_tmpdir: "/tmp"

Consequently, changing settings in `resources.default` will affect all
resource settings.

To modify resources for a rule, add the corresponding property in the
`resources` section under the rule name. For instance, to change
runtime, memory use, and threads for `gmap_map`, add

    resources:
      gmap_map:
        threads: 10
        runtime: 600
        mem_mb: 16000



## Testing

Test cases are in the subfolder `.test`. They are automatically
executed via continuous integration with [Github
Actions](https://github.com/features/actions). To run the tests, cd to
`.test` and issue

	snakemake --use-conda --conda-frontend mamba --show-failed-logs --cores 2 --conda-cleanup-pkgs cache -s ../workflow/Snakefile --wrapper-prefix file://$(pwd)/../workflow/wrappers

Once the test run has finished, create a report and view it:

	snakemake --cores 1 -s ../workflow/Snakefile --wrapper-prefix file://$(pwd)/../workflow/wrappers --report report.html
