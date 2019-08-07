Installation
============

Installation of Anaconda
------------------------

`MToolBox_snakemake`_, an update of the `MToolBox pipeline`_, is deployed in a conda environment, *i.e.* a virtual environment with all the needed tools/modules. Installing Anaconda is therefore essential, before installing the pipeline.

To this purpose, please follow instructions at http://docs.anaconda.com/anaconda/install/linux/ (hint: download the Anaconda installer in your personal directory with  `wget https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh`).

**Note on installation**: step 11 (verify installation by opening `anaconda-navigator`) is not compulsory. However, if you wish to do so, please make sure you have logged in the grid with either the `-X` or the `-Y` option.

Installation of MToolBox
------------------------

**Please note**: to `git clone` the repository you need to have an account on https://github.com/.

.. code-block:: bash
    
    # create directory and go there
    export pipelineDir="/path/to/MToolBox_snakemake"
    
    mkdir -p $pipelineDir
    cd $pipelineDir
    
    # update git
    conda install git
    
    # fetch repo
    git clone https://github.com/mitoNGS/MToolBox_snakemake.git
    
    # install environment
    cd MToolBox_snakemake
    conda env create \
    -n mtoolbox_snakemake \
    -f envs/environment.yaml

.. _`MToolBox_snakemake`: https://github.com/mitoNGS/MToolBox_snakemake 
.. _`MToolBox pipeline`: https://github.com/mitoNGS/MToolBox