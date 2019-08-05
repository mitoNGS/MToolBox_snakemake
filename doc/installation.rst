Installation
============

Installation of Anaconda
------------------------

The grid implementation of the MToolBox snakemake workflow, a side project of the MToolBox pipeline (https://github.com/mitoNGS/MToolBox), is deployed in a conda environment, *i.e.* a virtual environment with all the needed tools/modules. Installing Anaconda is therefore essential, before installing the pipeline.

To this purpose, please follow instructions at http://docs.anaconda.com/anaconda/install/linux/ (hint: download the Anaconda installer in your personal directory with  `wget https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh`).

**Note on installation**: step 11 (verify installation by opening `anaconda-navigator`) is not compulsory. However, if you wish to do so, please make sure you have logged in the grid with either the `-X` or the `-Y` option, *e.g.* `ssh -Y username@my-mgrid.mykopat.slu.se`.

Installation of the MToolBox-Ark workflow
-----------------------------------------

**Please note**: to `git clone` the repository you need to have an account on https://github.com/.

::
    
    # create directory and go there
    export pipelineDir="/nfs4/my-gridfront/mykopat-proj3/mykopat-hmtgen/heterobasidion_MToolBox"
    
    mkdir -p $pipelineDir
    cd $pipelineDir
    
    # update git
    conda install git
    
    # fetch repo
    git clone https://github.com/SLUBioinformaticsInfrastructure/MToolBox-Ark.git
    
    # install environment
    cd MToolBox-Ark
    conda env create \
    -n mtoolbox-ark \
    -f envs/environment.yaml
    
    # create folders needed by the workflow
    mkdir -p data/reads
    mkdir -p data/genomes
    mkdir -p logs/cluster_jobs
