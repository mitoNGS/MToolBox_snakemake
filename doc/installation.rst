Installation
============

Install Anaconda (or Miniconda) and mamba
-----------------------------------------

`MToolBox_snakemake`_, an update of the `MToolBox pipeline`_, is deployed in a conda environment, *i.e.* a virtual environment with all the needed tools/modules. Installing Anaconda (or Miniconda) is therefore essential, before installing the pipeline. *You may want to install Miniconda to keep your installation smaller* (since Anaconda comes with ~1500 prebuilt packages which you may not need!).

Please find instructions to install Anaconda `here <https://www.anaconda.com/products/distribution>`_ and Miniconda `here <https://docs.conda.io/en/latest/miniconda.html>`_.

**Please note**: MToolBox_snakemake has not been tested on Windows. In order to use MToolBox_snakemake on Windows, we strongly recommend to use the `Docker container <link>`.

Once Anaconda or Miniconda are installed, we also strongly recommend to install mamba to install the MToolBox_snakemake conda environment:

.. code-block:: bash
    
    conda install -c conda-forge mamba


Install MToolBox-snakemake
--------------------------

We recommend to download MToolBox-snakemake by cloning the official repo on `GitHub`_:

.. code-block:: bash
    
    # Pick a folder for your installation
    # and replace /path/to/MToolBox_snakemake with it
    cd /path/to/MToolBox_snakemake

    # fetch repo
    git clone https://github.com/mitoNGS/MToolBox_snakemake.git

Please note: you could also conveniently download MToolBox-snakemake at `this link`_, but by doing so you will miss the chance to easily integrate future updates!

Once you have cloned (or downloaded and unzipped) the repo, installing MToolBox should be as easy as running

.. code-block:: bash

   cd MToolBox_snakemake
   bash install.sh

The setup script ``install.sh`` will:

- install the ``mtoolbox`` conda environment with all the required dependencies
- create a command (``mtoolbox-activate``) which will be used to activate the MToolBox conda environment and add the folders of MToolBox executables and utilities to your ``PATH``.

.. note:: The ``install.sh`` will create the alias ``mtoolbox-activate`` in the ``~/.bash_profile`` file, so be sure it get sourced! 

.. _`MToolBox_snakemake`: https://github.com/mitoNGS/MToolBox_snakemake
.. _`MToolBox pipeline`: https://github.com/mitoNGS/MToolBox
.. _`GitHub`: https://github.com/
.. _`this link`: https://github.com/mitoNGS/MToolBox_snakemake/archive/master.zip
