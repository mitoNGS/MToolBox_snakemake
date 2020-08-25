Installation
============

Install Anaconda
----------------

`MToolBox_snakemake`_, an update of the `MToolBox pipeline`_, is deployed in a conda environment, *i.e.* a virtual environment with all the needed tools/modules. Installing Anaconda is therefore essential, before installing the pipeline.

To this purpose, please follow instructions at http://docs.anaconda.com/anaconda/install/linux/ (hint: download the Anaconda installer in your personal directory with  `wget https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh`).

Install MToolBox
----------------

MToolBox is hosted on `GitHub`_. You can get a copy of the repository by running this commands:

.. code-block:: bash

    cd /path/to/MToolBox_snakemake

    # update git
    conda install git

    # fetch repo
    git clone https://github.com/mitoNGS/MToolBox_snakemake.git

Installing MToolBox is as easy as running

.. code-block:: bash

   cd MToolBox_snakemake
   bash install.sh

The setup script ``install.sh`` will:

- install the ``mtoolbox`` conda environment with all the required dependencies
- create a command (``mtoolbox-activate``) which will be used to activate the MToolBox conda environment and add the folders of MToolBox executables and utilities to your ``PATH``.

.. _`MToolBox_snakemake`: https://github.com/mitoNGS/MToolBox_snakemake
.. _`MToolBox pipeline`: https://github.com/mitoNGS/MToolBox
.. _`GitHub`: https://github.com/
