assemblyeval-smk
================

The assemblyeval-smk_ workflow runs a suite of programs for evaluating
genome assembly quality. The workflow has rules to:

1. align transcripts to a reference sequence. Currently supported
   mappers are gmap_
2. given mapped transcripts in psl-format, estimate gene body coverage
   with the R package genecovr_
3. (WIP): run some of the steps for adding reads, coverage files,
   sequences to the blobtoolkit_ viewer

.. _assemblyeval-smk: https://github.com/percyfal/assemblyeval-smk
.. _genecovr: https://github.com/NBISweden/genecovr
.. _gmap: http://research-pub.gene.com/gmap/
.. _blobtoolkit: https://blobtoolkit.genomehubs.org/
