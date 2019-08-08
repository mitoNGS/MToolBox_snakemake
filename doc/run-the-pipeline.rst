Run MToolBox
============

MToolBox is made by several snakemake workflows which can be run independently. We provide wrappers for the most common tasks and analyses. Using these wrappers will save a lot of typing and headache for the lazy users (probably *you*). Cool, isn't it? :)
 
All the wrappers accepts snakemake **long** arguments (*e.g.* :code:`--dryrun` but not :code:`-n`) and parse automatically the `config.yaml` and `cluster.yaml` configuration files required by snakemake.

Setting up a working directory
------------------------------

Replace :code:`/path/to/MToolBox/dir/` with the MToolBox installation path and :code:`/path/to/analysis/dir` with the folder where you wish to run your analysis.

.. code-block:: bash
    
    export MTOOLBOX_DIR=/path/to/MToolBox/dir/
    
    cd /path/to/analysis/dir
    cp $MTOOLBOX_DIR/config.yaml .
    cp $MTOOLBOX_DIR/cluster.yaml .
    # create folders needed by the workflow
    mkdir -p data/reads
    mkdir -p data/genomes
    mkdir -p logs/cluster_jobs

Copy read datasets to `data/reads` and genome fasta files to `data/genomes`.

.. code-block:: bash
    
    # some code

How to run the MToolBox wrappers
--------------------------------

Running the wrappers is as simple as this:

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-<wrapper> <snakemake arguments>

*E.g.* if you want to run the MToolBox-variant-calling wrapper and print the commands it will execute, you can run

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling --printshellcmds

You can also display a graphical representation of the workflow by running

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling --printdag | dot -Tsvg | display

This will show the workflow in a browser. Alternatively, you can save the workflow representation in a file by running

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling --printdag | dot -Tsvg > workflow.svg

Available wrappers
------------------

- `MToolBox-variant-calling`_

MToolBox-variant-calling
^^^^^^^^^^^^^^^^^^^^^^^^

Performs QC, quality trimming of raw reads, read alignment, alignment filtering, variant calling. The final output is a VCF file.
