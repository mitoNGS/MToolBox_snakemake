Run MToolBox
============

MToolBox is made by several snakemake workflows which can be run independently. We provide wrappers for the most common tasks and analyses. Using these wrappers will save a lot of typing and headache for the lazy users (probably *you*). Cool, isn't it? :)
 
All the wrappers accepts snakemake **long** arguments and parse automatically the `config.yaml` and `cluster.yaml` configuration files required by snakemake.

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
