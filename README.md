# MToolBox on non-human genomes: Heterobasidion

<!-- TOC START min:1 max:4 link:true update:true -->
- [MToolBox on non-human genomes: Heterobasidion](#mtoolbox-on-non-human-genomes-heterobasidion)
	- [Installation](#installation)
		- [Installation of Anaconda](#installation-of-anaconda)
		- [Installation of the MToolBox workflow](#installation-of-the-mtoolbox-workflow)
		- [Copy/symlink data](#copysymlink-data)
	- [Running the pipeline](#running-the-pipeline)
		- [Activation of the conda environment](#activation-of-the-conda-environment)
		- [Compile configuration files](#compile-configuration-files)
			- [`data/analysis.tab`](#dataanalysistab)
			- [`data/reference_genomes.tab`](#datareference_genomestab)
		- [Run the whole workflow](#run-the-whole-workflow)
		- [Outputs](#outputs)
			- [Notes on outputs](#notes-on-outputs)
	- [Notes](#notes)
		- [Deactivation of the environment](#deactivation-of-the-environment)
	- [Graphical representation of the workflow](#graphical-representation-of-the-workflow)

<!-- TOC END -->


## Installation

### Installation of Anaconda

The grid implementation of the MToolBox snakemake workflow, a side project of the MToolBox pipeline (https://github.com/mitoNGS/MToolBox), is deployed in a conda environment, _i.e._ a virtual environment with all the needed tools/modules. Installing Anaconda is therefore essential, before installing the pipeline.

To this purpose, please follow instructions at http://docs.anaconda.com/anaconda/install/linux/ (hint: download the Anaconda installer in your personal directory with  `wget https://repo.continuum.io/archive/Anaconda3-2018.12-Linux-x86_64.sh`).

**Note on installation**: step 11 (verify installation by opening `anaconda-navigator`) is not compulsory. However, if you wish to do so, please make sure you have logged in the grid with either the `-X` or the `-Y` option, *e.g.* `ssh -Y username@my-mgrid.mykopat.slu.se`.

### Installation of the MToolBox workflow

Since the workflow and the input/output data will be in the same folder, it is strongly recommended to install the workflow in some subfolder of `/nfs4/my-gridfront/mykopat-proj3/mykopat-hmtgen/`. In this example this folder will be `/nfs4/my-gridfront/mykopat-proj3/mykopat-hmtgen/heterobasidion_MToolBox`.

**Please note**: to `git clone` the repository you need to have an account on https://github.com/.

```bash
# create directory and go there
export pipelineDir="/nfs4/my-gridfront/mykopat-proj3/mykopat-hmtgen/heterobasidion_MToolBox"

mkdir -p $pipelineDir
cd $pipelineDir

# update git
conda install git

# fetch repo
git clone https://github.com/domenico-simone/heterobasidion_mt.git

# install environment
cd heterobasidion_mt
conda env create \
-n heterobasidion_mt \
-f envs/environment.yaml

# create folders needed by the workflow
mkdir -p data/reads
mkdir -p data/genomes
mkdir -p logs
```

### Copy/symlink data

```bash
export projfolder="/nfs4/my-gridfront/mykopat-proj3/mykopat-hmtgen"

# reads
ln -s \
${projfolder}/hpar_raw_seq_data/*/*.fastq.gz \
data/reads

# nuclear genome
ln -s \
/nfs4/my-gridfront/mykopat-proj3/mykopat-hpadist/raw_sequence_data/reference/pb_121-1_polished_assembly.fasta \
data/genomes

# mt genomes
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

#### `data/analysis.tab`

Structure (strictly **tab-separated**):

| sample    | ref_genome_mt   | ref_genome_n   |
|:---------:|:---------------:|:--------------:|
| 87074_1 | 87074_1       | pb_121-1     |
| 87075_2 | 87074_1       | pb_121-1     |
| 87075_2 | 87075_2       | pb_121-1     |
| #A      | mt1           | n1           |

The workflow will execute as many analyses as the number of non-commented (i.e. not starting with `#`) lines in this file. Hint: instead of removing lines, you can comment them by prepending a `#` so that they will be skipped.
Every analysis is configured by one row in this table. Specifically:

- **sample** is the sample for which the read datasets will be used;
- **ref_genome_mt** is the reference mitochondrial genome used in the analysis;
- **ref_genome_n** is the reference mitochondrial genome used in the analysis.

*Eg*: The first row specifies that variant calling (with heteroplasmy) will be performed on sample 87074_1m using the mitochondrial reference genome 87074_1, by discarding those reads aligning on the nuclear reference genome pb_121-1.

#### `data/reference_genomes.tab`

Structure (strictly **tab-separated**):

| ref_genome_mt | ref_genome_n |      ref_genome_mt_file     |         ref_genome_n_file        |
|:-------------:|:------------:|:---------------------------:|:--------------------------------:|
|    87074_1    |   pb_121-1   | 87074_1-mtgenome-v270.fasta | pb_121-1_polished_assembly.fasta |
|    87075_2    |   pb_121-1   | 87075_2-mtgenome-v270.fasta | pb_121-1_polished_assembly.fasta |
|      mt1      |      n1      |          mt1.fasta          |             n1.fasta             |

This table contains explicit names for reference genome files used in the workflow. Names in the columns `ref_genome_mt` and `ref_genome_n` should be consistent with the ones in the same columns in the `data/analysis.tab` table. Genome files must be located in the `data/genomes` folder.

### Run the whole workflow

```bash
nohup \
snakemake -rp \
-j 100 --cluster-config cluster.yaml \
--cluster 'qsub -V -l h_rt={cluster.time} -l h_vmem={cluster.vmem} -pe smp {cluster.threads} -cwd -j y -o {cluster.stdout}' &> logs/nohup_MToolBox.log &
```

The file `logs/nohup_MToolBox.log` contains info about the whole run. The `logs` folder will also contain logs for each step of the workflow.

### Outputs

Folders created during the workflow execution: `results`, `gmap_db`, `logs`.

- `results` folder tree

```
results/
├── OUT_87074_1_87074_1_pb_121-1
│   ├── map
│   │   ├── OUT.bam
│   │   ├── outmt1.fastq.gz
│   │   ├── outmt2.fastq.gz
│   │   ├── outmt.fastq.gz
│   │   ├── outmt.sam.gz
│   │   ├── outP.sam.gz
│   │   ├── OUT.sam.gz
│   │   ├── OUT-sorted.bam
│   │   └── outS.sam.gz
│   ├── variant_calling
│   │   ├── OUT-mt_table.txt
│   │   └── OUT-sorted.pileup
│   └── vcf.vcf
└── OUT_87075_2_87075_2_pb_121-1
    ├── map
    │   ├── OUT.bam
    │   ├── outmt1.fastq.gz
    │   ├── outmt2.fastq.gz
    │   ├── outmt.fastq.gz
    │   ├── outmt.sam.gz
    │   ├── outP.sam.gz
    │   ├── OUT.sam.gz
    │   ├── OUT-sorted.bam
    │   └── outS.sam.gz
    ├── variant_calling
    │   ├── OUT-mt_table.txt
    │   └── OUT-sorted.pileup
    └── vcf.vcf
```

#### Notes on outputs

For each sample:

- `OUT-sorted.bam` contains filtered alignments (in binary format) used for variant calling, useful for displaying alignment data (reads overlapping a specific genome region, coverage depth at specific sites etc) in external tools, _e.g._ IGV.
- `OUT-sorted.pileup` contains a text-format report of SNP calling. It can suggest the presence of indels, but it does not report them. More details at https://en.wikipedia.org/wiki/Pileup_format.
- `vcf.vcf` contains genotype information for each sample. The genotype is reported following the structure reported in the FORMAT field GT:DP:HF:CILOW:CIUP, where each subfield is explained as follows:

|  Subfield |                Description                |
|:---------:|:-----------------------------------------:|
|   **GT**  |                  genotype                 |
|   **DP**  |      reads covering the REF position      |
|   **HF**  |   heteroplasmy fraction of ALT allele(s)  |
| **CILOW** | lower limit of the HF confidence interval |
|  **CIUP** | upper limit of the HF confidence interval |

For example, if sample A shows this genotype information:  

| #CHROM  | POS   | ID | REF | ALT     | QUAL | FILTER | INFO        | FORMAT              | 87074_1                                     |
|---------|-------|----|-----|---------|------|--------|-------------|---------------------|---------------------------------------------|
| Contig1 | 26169 | .  | TT  | T       | .    | PASS   | AC=1;AN=2   | GT:DP:HF:CILOW:CIUP | 0/1:6895:0.01:0.008:0.013                   |
| Contig1 | 26158 | .  | AT  | ATTTT,A | .    | PASS   | AC=1,1;AN=3 | GT:DP:HF:CILOW:CIUP | 0/1/2:5553:0.001,0.016:0.0,0.014:0.002,0.02 |

- at site 26169, sample A has a 1T deletion with a heteroplasmy frequency of 0.01 (1%). That site is covered by 6895 reads. The heteroplasmy frequency is calculated with a confidence interval. In this case, the confidence interval ranges from 0.008 and 0.013 (0.8-0.13%).

- at site 26158, sample A has two variants: an insertion of three Ts (with HF = 0.001) and a 1T deletion with HF of 0.016. That site is covered by 5553 reads. The HF confidence intervals for each variant are 0-0.002 and 0.014-0.02, respectively.

Please find more info about the VCF format here: http://www.internationalgenome.org/wiki/Analysis/vcf4.0/ and more info about the calculation of HF and the related confidence interval in the original MToolBox paper: https://www.ncbi.nlm.nih.gov/pubmed/25028726.

## Notes

### Deactivation of the environment

```bash
# To deactivate an active environment, use
conda deactivate
```

- if the same cluster job is runned again, its log will be appended to the existing file [need to change this behaviour]

## Graphical representation of the workflow

A graphical representation of the workflow can be obtained by running

```bash
snakemake --dag | dot -Tsvg > my_workflow.svg
```

The graph in file `my_workflow.svg` will report all the workflow steps (for each sample in the `analysis.tab` configuration file). Steps in dashed lines are to be run (because their outputs are not present), whereas outputs for steps in solid lines are already present. A graphical representation of the workflow as per the `analysis.tab` file in this repo is reported as follows.

![workflow](workflow.png)
