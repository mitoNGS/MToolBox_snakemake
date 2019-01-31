# MToolBox on non-human genomes: Heterobasidion

<!-- TOC START min:1 max:3 link:true update:true -->
- [MToolBox on non-human genomes: Heterobasidion](#mtoolbox-on-non-human-genomes-heterobasidion)
	- [Installation](#installation)
		- [Installation of Anaconda](#installation-of-anaconda)
		- [Installation of the MToolBox workflow](#installation-of-the-mtoolbox-workflow)
		- [Copy/symlink data](#copysymlink-data)
	- [Running the pipeline](#running-the-pipeline)
		- [Activation of the conda environment](#activation-of-the-conda-environment)
		- [Compile configuration files](#compile-configuration-files)
		- [Run the whole pipeline](#run-the-whole-pipeline)
	- [Notes](#notes)
		- [Deactivation of the environment](#deactivation-of-the-environment)
	- [Graphical representation of the workflow](#graphical-representation-of-the-workflow)

<!-- TOC END -->


## Installation

### Installation of Anaconda

The grid implementation of the MToolBox pipeline is deployed in a conda environment, _i.e._ a virtual environment with all the needed tools/modules. Installing conda is therefore essential, before installing the pipeline.

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

# create folders needed by the workflow
mkdir -p data/reads
mkdir -p data/genomes
```

### Copy/symlink data

```bash
export projfolder="/nfs4/my-gridfront/mykopat-proj3/mykopat-hmtgen"

# reads
ln -s \
${projfolder}/hpar_raw_seq_data/*/*.fastq.gz \
data/reads

# genomes
ln -s \
${projfolder}/NOVOPlasty_assemblys/*/*.fasta \
data/genomes
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

### Compile configuration files

- `data/analysis.tab`
- `data/reference_genomes.tab`

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

## Graphical representation of the workflow

![workflow](workflow.png)
