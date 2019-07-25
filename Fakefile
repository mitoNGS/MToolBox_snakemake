import pandas as pd
import os, re, sys, time, gzip, bz2, subprocess
from Bio import SeqIO
import resource
import numpy as np
#import sqlite3
from sqlalchemy import create_engine
from modules.mtVariantCaller import *
from modules.BEDoutput import *

# def fastqc_filtered_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_filtered", ext=".fastq.gz"):
#     fastqc_out = []
#     for s in analysis_tab["sample"]:
#         # fastqc_html_1 and fastqc_html_2 could be strings, but
#         # keep them as lists in case one sample has multiple datasets
#         #L = set(os.path.join(infolder))
#         #glob.glob("{in}/{}*R1{}.format()")
#         # expand("results/fastqc_filtered/{sample}")
#         # fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
#         # fastqc_out.extend(fastq_html_1)
#         #######
#         for read_type in ["R1", "R2", "U"]:
#             fastqc_out.append(os.path.join(outfolder, "{sample}.{read_type}_fastqc.html".format(sample=s, read_type=read_type)))
#     return fastqc_out

def fastqc_raw_outputs(datasets_tab, infolder="data/reads", outfolder="results/fastqc_raw", ext=".fastq.gz"):
    fastqc_out = []
    for i,l in datasets_tab.iterrows():
        fastqc_out.append(os.path.join(outfolder, l["R1"].replace(ext, "_fastqc.html")))
        fastqc_out.append(os.path.join(outfolder, l["R2"].replace(ext, "_fastqc.html")))
    return fastqc_out

datasets_tab = pd.read_table("data/datasets.tab", sep = "\t", comment='#')
analysis_tab = pd.read_table("data/analysis.tab", sep = "\t", comment='#')

wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in list(set(analysis_tab['sample']))]),
    ref_genome_mt = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_mt']))]),
    ref_genome_n = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_n']))])

#print(fastqc_raw_outputs(datasets_tab))

rule all:
    input: fastqc_raw_outputs(datasets_tab)

rule fastqc_raw:
    input:
        R1 = "data/reads/{sample}_{adapter}_{lane}_R1_001.fastq.gz",
        R2 = "data/reads/{sample}_{adapter}_{lane}_R2_001.fastq.gz"
        # R1 = lambda wildcards: read_datasets_inputs(sample="{sample}".format(sample=wildcards.sample), read_type="1", input_folder="data/reads"),
        # R2 = lambda wildcards: read_datasets_inputs(sample="{sample}".format(sample=wildcards.sample), read_type="2", input_folder="data/reads"),
    output:
        html_report_R1 = "results/fastqc_raw/{sample}_{adapter}_{lane}_R1_001_fastqc.html",
        html_report_R2 = "results/fastqc_raw/{sample}_{adapter}_{lane}_R2_001_fastqc.html",
        #html_report_R1 = "results/fastqc_raw/{sample}_R1_fastqc.html",
        #html_report_R2 = "results/fastqc_raw/{sample}_R2_fastqc.html",
        # html_report_R1 = lambda wildcards: "results/fastqc_raw/{outR1}".format(outR1 = os.path.split({input.R1})[1].replace(".fastq.gz", "_fastqc.html")),
        # html_report_R2 = lambda wildcards: "results/fastqc_raw/{outR2}".format(outR2 = os.path.split({input.R2})[1].replace(".fastq.gz", "_fastqc.html"))
        #logFile = os.path.join(config["proj_dirs"]["logs"], "fastqc_raw.log")
    params:
        outDir = "results/fastqc_raw/",
    threads:
        2
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of raw read files {input} with {version}, {wildcards}"
    # log:
    #     "logs/fastqc_raw/{sample}.log"
    shell:
        """
        cd $(dirname {input.R1})
        mkdir -p {wildcards.sample}
        ln -sf `pwd`/$(basename {input.R1}) {wildcards.sample}/{wildcards.sample}_R1.fastq.gz
        ln -sf `pwd`/$(basename {input.R2}) {wildcards.sample}/{wildcards.sample}_R2.fastq.gz
        cd -
        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} data/reads/{wildcards.sample}/{wildcards.sample}_R1.fastq.gz data/reads/{wildcards.sample}/{wildcards.sample}_R2.fastq.gz > {log}
        rm -R data/reads/{wildcards.sample}
        """

rule merge_fastqc:
    input: