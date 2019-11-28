import bz2
import gzip
import os
import re
import resource
import shutil
import subprocess
import sys
import time

from Bio import SeqIO
import numpy as np
import pandas as pd
from sqlalchemy import create_engine

# TODO: Not sure this is the right way to deal with this
for path in sys.path:
    if "snakefiles" in path:
        sys.path.append(path.replace("/snakefiles", ""))

from modules.config_parsers import get_haplo_prediction_files
from modules.mt_classifier import main_mt_hpred, write_output


source_dir = "/".join(os.path.dirname(workflow.snakefile).split("/")[:-1])

analysis_tab = pd.read_table("data/analysis.tab", sep = "\t", comment='#')
reference_tab = (pd.read_table("data/reference_genomes.tab", sep = "\t", comment='#')
                 .set_index("ref_genome_mt", drop=False))
datasets_tab = pd.read_table("data/datasets.tab", sep = "\t", comment='#')

include: "variant_calling.snakefile"

configfile: "config.yaml"
res_dir = config["results"]
log_dir = config["log_dir"]

wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in list(set(analysis_tab['sample']))]),
    ref_genome_mt = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_mt']))]),
    ref_genome_n = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_n']))])

rule all_haplo_prediction:
    input:
        get_haplo_prediction_files(analysis_tab)

rule main_mt_hpred:
    input:
        single_fasta = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        haplo_pred = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.csv"
    params:
        basename = lambda wildcards: "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}".format(
            sample=wildcards.sample, ref_genome_mt=wildcards.ref_genome_mt,
            ref_genome_n=wildcards.ref_genome_n
        )
    run:
        sc, contig_seq_diff, contig_mhcs_seq_diff, contig_rcrs_seq_diff, mergedtables = main_mt_hpred(
            contig_file=contig_file, muscle_exe=muscle_exe, basename=basename,
            best_results_file=best_results_file, data_file=data_file
        )
        write_output(sc, contig_seq_diff.diff_list, contig_mhcs_seq_diff.diff_list,
                     contig_rcrs_seq_diff.diff_list, mergedtables, basename)

