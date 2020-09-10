Run MToolBox
============

MToolBox is made by several snakemake workflows which can be run independently. We provide wrappers for the most common tasks and analyses. Using these wrappers will save a lot of typing and headache for the lazy users (probably *you*). Cool, isn't it? :)
 
All the wrappers accepts snakemake arguments and parse automatically the `config.yaml` configuration file required by snakemake.

An overview
-----------

You are going to analyse one or more **samples**, be represented by one or more read **datasets** (*ie*, libraries). 
For this purpose, you're going to provide a **reference mitochondrial genome**. Choose it carefully, as your final results will be based on it! You also have to pick a **reference nuclear genome** that will be used as reference to filter out ambiguous reads.

If you are interested in performing functional annotation of mt variants, you will also explicitly provide **species**.

Now you'll be wondering: *how do I tell all these things to the pipeline?* In the following sections, we'll see how to do that through a handful of configuration files!

Setting up a working directory
------------------------------

**Note:** Replace :code:`/path/to/MToolBox/dir/` with the MToolBox installation path and :code:`/path/to/analysis/dir` with the folder where you wish to run your analysis.

.. code-block:: bash
    
    export MTOOLBOX_DIR=/path/to/MToolBox/dir/
    
    cd /path/to/analysis/dir
    
    # create folders needed by the workflow
    mkdir -p data/reads
    mkdir -p data/genomes
    mkdir -p logs/cluster_jobs
    
    # copy configuration files you will edit later
    cp $MTOOLBOX_DIR/config.yaml .
    cp $MTOOLBOX_DIR/cluster.yaml .
    cp $MTOOLBOX_DIR/data/*.tab data

At this point, if you run the command :code:`tree` the structure of your directory should look like

.. code-block:: bash
    
    .
    ├── cluster.yaml
    ├── config.yaml
    ├── data
    │   ├── analysis.tab
    │   ├── datasets.tab
    │   ├── genomes
    │   ├── reads
    │   └── reference_genomes.tab
    └── logs
        └── cluster_jobs

Compiling configuration files
-----------------------------

An MToolBox-snakemake run is managed with these configuration files: 

- data/analysis.tab
- data/reference_genomes.tab
- data/datasets.tab
- config.yaml
- cluster.yaml

*Sounds a pain, huh?* The good news is that they will help you in setting up and keeping track of your analyses very efficiently. Plus, the config.yaml and cluster.yaml files should work the way they are, with little to no edit needed. **Please read the "Notes on configuration files" at the end of this section**.

Let's see how to compile the configuration files in detail.

- **data/analysis.tab**

For each sample you are going to analyse, in this table you provide info about which mitochondrial and nuclear reference genomes to use. Example:

+----------+---------------+-----------------+
| sample   | ref_genome_mt | ref_genome_n    |
+==========+===============+=================+
| sample_1 | NC_001323.1   | GCF_000002315.5 |
+----------+---------------+-----------------+
| sample_2 | NC_001323.1   | GCF_000002315.5 |
+----------+---------------+-----------------+

In this example, the first row specifies that variant calling will be performed on **sample_1** using the mitochondrial reference genome **NC_001323.1**, by discarding those reads aligning on the nuclear reference genome **GCF_000002315.5**. Please note that the names used in this table will be used in the workflow execution and are case-sensitive. Actual files related to samples and reference genomes will be provided in the **data/reference_genomes.tab** and in the **data/datasets.tab** files.

- **data/reference_genomes.tab**

Structure (strictly **tab-separated**):

+---------------+-----------------+--------------------+-----------------------+---------+
| ref_genome_mt | ref_genome_n    | ref_genome_mt_file | ref_genome_n_file     | species |
+===============+=================+====================+=======================+=========+
| NC_001323.1   | GCF_000002315.5 | NC_001323.1.fasta  | GCF_000002315.5.fasta | ggallus |
+---------------+-----------------+--------------------+-----------------------+---------+

This table contains explicit names for reference genome files used in the workflow. Names in the columns **ref_genome_mt** and **ref_genome_n** must be consistent with the ones in the same columns in the **data/analysis.tab** table. **Genome files must be located in the data/genomes folder**.

The name in the column **species** should be one of the `species available in mtoolnote`_ for variant functional annotation. 

- **data/datasets.tab**

Fill this table with as many read (paired) datasets are available per sample. Each read dataset will be processed independently and merged with the others from the same sample before the variant calling stage. **Read dataset files must be located in the data/reads folder**.

Example:

+----------+---------+--------------------------+--------------------------+
| sample   | library | R1                       | R2                       |
+==========+=========+==========================+==========================+
| sample_1 | 1       | sample_1_R1_001.fastq.gz | sample_1_R2_001.fastq.gz |
+----------+---------+--------------------------+--------------------------+
| sample_1 | 2       | sample_1_R1_002.fastq.gz | sample_1_R2_002.fastq.gz |
+----------+---------+--------------------------+--------------------------+
| sample_2 | 1       | sample_2_R1.fastq.gz     | sample_2_R2.fastq.gz     |
+----------+---------+--------------------------+--------------------------+

In this case, sample_1 is represented by two PE libraries, while sample_2 is represented by one.

- **config.yaml**

This file contains basic configuration for the whole workflow. Default configuration should fit most cases; you might want to check the `mark_duplicates` option (which removes duplicate reads with Picard MarkDuplicates) and set it to True or False, depending on your needs.
TODO: add realign indels option

- **cluster.yaml**

TODO: add stuff

A recap
^^^^^^^

.. figure:: img/MToolBox_conf_files.png
    :align: center
    :alt: alternate text
    :figclass: align-center

    An overview of MToolBox-snakemake configuration files

How to run the MToolBox wrappers
--------------------------------

Running the wrappers is as simple as this:

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-<wrapper> <snakemake arguments>

*E.g.* if you want to run the MToolBox-variant-calling wrapper and print the commands it will execute, you can run

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling -p

You can also display a graphical representation of the workflow by running

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling --dag | display

This will show the workflow in a browser. Alternatively, you can save the workflow representation in a file by running

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling --dag > workflow.svg

Available wrappers
------------------

- `MToolBox-variant-calling`_

MToolBox-variant-calling
^^^^^^^^^^^^^^^^^^^^^^^^

Performs QC, quality trimming of raw reads, read alignment, alignment filtering, variant calling. The final output is a VCF file.

.. _`species available in mtoolnote`: https://github.com/mitoNGS/mtoolnote#features