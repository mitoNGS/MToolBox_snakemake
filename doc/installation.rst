Installation
============

Install Anaconda
----------------

`MToolBox_snakemake`_, an update of the `MToolBox pipeline`_, is deployed in a conda environment, *i.e.* a virtual environment with all the needed tools/modules. Installing Anaconda is therefore essential, before installing the pipeline.

To this purpose, please follow instructions at http://docs.anaconda.com/anaconda/install/linux/ (hint: download the Anaconda installer in your personal directory with  `wget https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh`).

**Note on installation**: step 11 (verify installation by opening `anaconda-navigator`) is not compulsory. However, if you wish to do so, please make sure you have logged in the grid with either the `-X` or the `-Y` option.

Install MToolBox
----------------

MToolBox is hosted on `GitHub`_. You can get a copy by running this commands:

.. code-block:: bash

    cd /path/to/MToolBox_snakemake

    # update git
    conda install git

    # fetch repo
    git clone https://github.com/mitoNGS/MToolBox_snakemake.git

The MToolBox repo comes with a setup script `install.sh`, which will:

- install the `mtoolbox` conda environment
- install the `bamUtils` suite (a third-party tool used in one of the steps of the pipeline)
- create a command (`mtoolbox_activate`) which will be used to activate the `mtoolbox` conda environment and add the folders of MToolBox executables and utilities to your `PATH`.

.. code-block:: bash

    cd MToolBox_snakemake
    bash install.sh

.. _`MToolBox_snakemake`: https://github.com/mitoNGS/MToolBox_snakemake
.. _`MToolBox pipeline`: https://github.com/mitoNGS/MToolBox
.. _`GitHub`: https://github.com/
