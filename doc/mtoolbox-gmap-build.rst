.. _mtoolbox_genomedb_build:

Get your reference genome(s)
============================

This workflow is aimed at getting your favourite reference genome(s) and build a reference genome database for GMAP/GSNAP that you can re-use for subsequent analyses. We provide pre-built reference genome dbs will be downloaded for species whose annotation is available through mtoolnote.

How it works
------------

A standard command line for the MToolBox-genomedb-build workflow is

.. code-block:: bash
    
    MToolBox-genomedb-build --config ref_organism=ggallus,human

which will create GMAP reference databases for ggallus and human. If the provided ref_organisms are available as mtoolnote species, you will just download the GMAP reference db from our magical horn of plenty (at Zenodo). Otherwise, please go ahead and read the next section!

Build a reference genome db for non-mtoolnote species
-----------------------------------------------------

If the ref_organism you provide is not an mtoolnote species, you have to carefully compile the :code:`reference_genomes.tab` file. The structure of the file is

+---------------+-----------------+--------------------+---------------------------------+--------------------------+--------------------------+
| ref_genome_mt | ref_genome_n    | ref_genome_mt_file | ref_genome_n_file               | species                  | ref_organism             |
+===============+=================+====================+=================================+==========================+==========================+
| NC_001224.1   | GCF_000146045.2 | NC_001224.1.fa     | GCF_000146045.2_R64_genomic.fna | scerevisiae              | scerevisiae_2            |
+---------------+-----------------+--------------------+---------------------------------+--------------------------+--------------------------+
| NC_002333.2   | GCF_000000175.5 |                    |                                 | drerio                   | drerio_2                 |
+---------------+-----------------+--------------------+---------------------------------+--------------------------+--------------------------+
| NC_008282.1   | GCA_000184455.3 |                    |                                 | aspergillus_oryzae_RIB40 | aspergillus_oryzae_RIB40 |
+---------------+-----------------+--------------------+---------------------------------+--------------------------+--------------------------+

