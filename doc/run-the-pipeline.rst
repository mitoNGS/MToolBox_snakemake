Run the pipeline
================

MToolBox is made by several snakemake workflows which can be run independently. We provide wrappers for the most common tasks and analyses, which automatically parse the `config.yaml` and `cluster.yaml` configuration files required by snakemake. Using this wrappers will save a lot of typing and headache for the lazy users (probably *you*). Cool, isn't it? :)
 
All the wrappers accepts snakemake *long* arguments.

Available wrappers
------------------

- `MToolBox-variant-calling`_

MToolBox-variant-calling
^^^^^^^^^^^^^^^^^^^^^^^^

Performs QC, quality trimming of raw reads, read alignment, alignment filtering, variant calling. The final output is a VCF file.
