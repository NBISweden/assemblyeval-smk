assemblyeval-smk
================

The assemblyeval-smk_ workflow runs a suite of programs for evaluating
genome assembly quality. This analysis is based on commit version {{
snakemake.config["__workflow_commit__"] }}.

The analysis can be rerun with the following command:

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --use-conda --wrapper-prefix {{ "file://" + snakemake.config["__workflow_basedir__"] + "/wrappers" }}
{% else %}
   snakemake -j 1 --use-conda -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile  --wrapper-prefix {{ "file://" + snakemake.config["__workflow_basedir__"] + "/wrappers" }}
{% endif %}

and the report

.. code-block:: shell

   cd {{ snakemake.config["__workflow_workdir__"] }}
{% if snakemake.config["__workflow_basedir__"] == snakemake.config["__workflow_workdir__"] %}
   snakemake -j 1 --wrapper-prefix {{ "file://" + snakemake.config["__workflow_basedir__"] + "/wrappers" }} --report report.html
{% else %}
   snakemake -j 1 -s {{ snakemake.config["__workflow_basedir__"] }}/Snakefile --wrapper-prefix {{ "file://" + snakemake.config["__workflow_basedir__"] + "/wrappers" }}  --report report.html
{% endif %}

Note the use of the option `--wrapper-prefix`. All rules make use of
`Snakemake wrappers`_, where some are based on existing wrappers in the
snakemake wrapper repository, whereas others are custom-made wrappers
for this workflow. Briefly, a wrapper consists of a metadata file, an
environment file that defines software dependencies, and a python
wrapper script that determines how the application is run. See {{
"file://" + snakemake.config["__workflow_basedir__"] +
"/wrappers/bio/gmap/map" }} for an example.


Workflow summary
-----------------

The workflow has rules to:

- align transcripts to a reference sequence. Currently supported
  mappers are gmap_
- given mapped transcripts in psl-format, estimate gene body coverage
  with the R package genecovr_
- run quast_, jellyfish_, and busco_
- summarize quality metrics with MultiQC_
- produce kmer-plots of assembly and reads
- (WIP): run some of the steps for adding reads, coverage files,
  sequences to the blobtoolkit_ viewer



Data organization
=================

.. code-block:: text

   {{ snakemake.config["__workflow_workdir__"] }}/                                <- top-level project folder
   |
   ├── config                   <- configuration directory for Snakemake and other things
   │
   ├── data
   │   ├── external             <- data from third party sources
   │   ├── interim              <- Intermediate data that can be safely deleted
   │   ├── metadata             <- metadata describing raw data files
   │   ├── processed            <- Final processed data used for analyses
   │   └── raw                  <- The original immutable data dump to be treated as read-only.
   │
   ├── logs                     <- Collection of log outputs, e.g. from cluster managers
   │
   ├── reports                  <- Generated analyses and articles as html, pdf and more, including multiqc.html
   │   └── figures              <- Graphics for use in reports.
   │
   └── results                  <- Final results for sharing with collaborators, typically derived



MultiQC report
=================

The `MultiQC report`_ collects results from quast_, jellyfish_, busco_,
and other QC programs for which `MultiQC tools`_ exist.

.. _assemblyeval-smk: https://github.com/percyfal/assemblyeval-smk
.. _genecovr: https://github.com/NBISweden/genecovr
.. _gmap: http://research-pub.gene.com/gmap/
.. _blobtoolkit: https://blobtoolkit.genomehubs.org/
.. _MultiQC: https://multiqc.info/
.. _quast: http://bioinf.spbau.ru/quast
.. _jellyfish: http://www.genome.umd.edu/jellyfish.html
.. _busco: https://busco.ezlab.org/
.. _Snakemake wrappers: https://snakemake-wrappers.readthedocs.io/en/stable/
.. _MultiQC report: ./reports/multiqc.html
.. _MultiQC tools: https://multiqc.info/#supported-tools
