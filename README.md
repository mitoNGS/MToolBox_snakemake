# MToolBox on non-human genomes: Heterobasidion

## Installation

### Installation of Anaconda 

The grid implementation of the Mykopat pipeline is deployed in a conda environment, i.e. a virtual environment with all the needed tools/modules (including hmmer and clustal-omega). Installing conda is therefore essential, before installing the pipeline.

Follow instructions at http://docs.anaconda.com/anaconda/install/linux/

### Installation of the MToolBox workflow

```bash
# update git
conda install git

# fetch repo
git clone https://github.com/domenico-simone/heterobasidion_mt.git

# install environment
conda env create \
-n heterobasidion_mt \
-f envs/environment.yaml
```

## Running the pipeline

### Activation of the conda environment

```bash
# To activate this environment, use
source activate heterobasidion_mt

# load appropriate modules (their conda packages have problems)
module load samtools/1.8
module load gsnap
```

### Run the whole pipeline

```bash
nohup \
snakemake -rp \
-j 100 --cluster-config cluster.yaml \
--cluster 'qsub -V -l h_rt={cluster.time} -l h_vmem={cluster.vmem} -pe smp {cluster.threads} -cwd -j y -o {cluster.stdout}' &> logs/nohup_MToolBox.log &
```

The file `logs/nohup_MToolBox.log` contains info about the whole run. The `logs` folder will also contain logs for each step of the workflow.

## Notes

### Deactivation of the environment

```bash
# To deactivate an active environment, use
conda deactivate
```

- if the same cluster job is runned again, its log will be appended to the existing file [need to change this behaviour]