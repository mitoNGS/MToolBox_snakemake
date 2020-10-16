import bz2
import gzip
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time
import yaml

from Bio import SeqIO
import numpy as np
import pandas as pd

current_dir = os.getcwd()
for path in sys.path:
    if "snakefiles" in path:
        sys.path.append(path.replace("/snakefiles", ""))
        rootdir = path.replace("/snakefiles", "")
        snakedir = path

from modules.BEDoutput import bed_output, fasta_output
from modules.config_parsers import (
    fastqc_outputs, get_bed_files, get_datasets_for_symlinks,
    get_fasta_files, get_genome_files, get_genome_single_vcf_files,
    get_genome_single_vcf_index_files, get_genome_vcf_files, get_inputs_for_rule_map_nuclear_MT_SE,
    get_mt_genomes, get_mt_fasta, get_sample_bamfiles, get_symlinks, parse_config_tabs, parse_config_tab
)
from modules.filter_alignments import filter_alignments, cat_alignments
from modules.general import (
    check_tmp_dir, gapped_fasta2contigs, get_seq_name, sam_to_fastq, sam_cov_handle2gapped_fasta,
    trimmomatic_input, sam_to_ids, run_seqtk_subset, get_mt_AC, get_SAM_header
)
from modules.genome_db import (
    check_ref_organism, make_ref_organism_dict
)
from modules.mtVariantCaller import mtvcf_main_analysis, VCFoutput

source_dir = Path(os.path.dirname(workflow.snakefile)).parent
localrules: merge_VCF, index_VCF, symlink_libraries, symlink_libraries_uncompressed

genome_db_data_file = os.path.join(rootdir, "data/genome_dbs.yaml")
reference_tab_file =  os.path.join(rootdir, "data/reference_genomes.tab")
with open(genome_db_data_file) as file:
    genome_db_data = yaml.full_load(file)

# keep in mind there's a parse_config function that calls these three functions
# analysis_tab, reference_tab, datasets_tab = parse_config_tabs(analysis_tab_file="data/analysis.tab", reference_tab_file="data/reference_genomes.tab", datasets_tab_file="data/datasets.tab")
analysis_tab = parse_config_tab(tab_file="data/analysis.tab", index=["sample"])
reference_tab = parse_config_tab(tab_file=reference_tab_file, index=["ref_organism"])
datasets_tab = parse_config_tab(tab_file="data/datasets.tab", index=["sample", "library"])

configfile: "config.yaml"
res_dir = config["results"]
map_dir = config["map_dir"]
log_dir = config["log_dir"]
gmap_db_dir = os.path.join(rootdir, config["map"]["gmap_db_dir"])
keep_orphans = config["keep_orphans"]
genome_fasta_dir = os.path.join(rootdir, "data/genomes")

ref_organism_config, analysis_tab = check_ref_organism(config=config, analysis_tab=analysis_tab)

ref_organism_dict = make_ref_organism_dict(ref_organism_config=ref_organism_config,
                                            gmap_db_dir=gmap_db_dir, genome_db_data=genome_db_data,
                                            genome_fasta_dir=genome_fasta_dir, reference_tab=reference_tab)

wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in list(set(analysis_tab['sample']))]),
    ref_organism = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_organism']))]),
    # ref_genome_mt = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_mt']))]),
    # ref_genome_n = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_n']))])

subworkflow genome_db:
    workdir:
        #snakedir
        current_dir
    snakefile:
        os.path.join(rootdir, "snakefiles", "genome_db.snakefile")
    configfile:
        os.path.join(current_dir, "config.yaml")

rule all:
    input:
        get_symlinks(datasets_tab, analysis_tab=analysis_tab, infolder="data/reads",
                     outfolder="data/reads"),
        fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="raw"),
        fastqc_outputs(datasets_tab, analysis_tab=analysis_tab, out="filtered"),
        get_genome_vcf_files(analysis_tab),
        get_bed_files(analysis_tab),
        get_fasta_files(analysis_tab)

rule symlink_libraries:
    input:
        R1 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         library=wildcards.library, d="R1"),
        R2 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         library=wildcards.library, d="R2")
    output:
        R1 = "data/reads/{sample}_{library}.R1.fastq.gz",
        R2 = "data/reads/{sample}_{library}.R2.fastq.gz",
    shell:
        """
        cd data/reads/
        ln -sf $(basename {input.R1}) $(basename {output.R1})
        ln -sf $(basename {input.R2}) $(basename {output.R2})
        """

rule symlink_libraries_uncompressed:
    input:
        R1 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         library=wildcards.library, d="R1"),
        R2 = lambda wildcards: get_datasets_for_symlinks(datasets_tab, sample=wildcards.sample,
                                                         library=wildcards.library, d="R2")
    output:
        R1 = "data/reads/{sample}_{library}.R1.fastq",
        R2 = "data/reads/{sample}_{library}.R2.fastq",
    shell:
        """
        cd data/reads/
        ln -sf $(basename {input.R1}) $(basename {output.R1})
        ln -sf $(basename {input.R2}) $(basename {output.R2})
        """

rule fastqc_raw:
    input:
        R1 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library)[0],
        R2 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library)[1]
    output:
        html_report_R1 = "results/fastqc_raw/{sample}_{library}.R1_fastqc.html",
        html_report_R2 = "results/fastqc_raw/{sample}_{library}.R2_fastqc.html",
    params:
        outDir = "results/fastqc_raw/",
    threads:
        2
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of raw read files {input} with {version}, {wildcards}"
    log:
        "logs/fastqc_raw/{sample}_{library}.log"
    #conda: "envs/environment.yaml"
    shell:
        """
        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} {input} &> {log}
        """

rule fastqc_filtered:
    input:
        out1P = "data/reads_filtered/{sample}_{library}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{library}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{library}_qc_U.fastq.gz",
    output:
        html_report_R1 = "results/fastqc_filtered/{sample}_{library}_qc_R1_fastqc.html",
        html_report_R2 = "results/fastqc_filtered/{sample}_{library}_qc_R2_fastqc.html",
        html_report_U = "results/fastqc_filtered/{sample}_{library}_qc_U_fastqc.html",
    params:
        outDir = "results/fastqc_filtered/"
    threads:
        3
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of filtered read files {input} with {version}"
    log:
        "logs/fastqc_filtered/{sample}_{library}.log"
    #conda: "envs/environment.yaml"
    shell:
        """

        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} {input} &> {log}

        """

rule trimmomatic:
    """ QCing and cleaning reads """
    params:
        java_cmd = config['read_processing']['trimmomatic']['java_cmd'],
        #jar_file = config['read_processing']['trimmomatic']['jar_file'],
        mem = config['read_processing']['trimmomatic']['java_vm_mem'],
        options = config['read_processing']['trimmomatic']['options'],
        processing_options = config['read_processing']['trimmomatic']['processing_options'],
        out1P = "data/reads_filtered/{sample}_{library}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{library}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{library}_qc_1U.fastq.gz",
        out2U = "data/reads_filtered/{sample}_{library}_qc_2U.fastq.gz"
    input:
        R1 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library)[0],
        R2 = lambda wildcards: trimmomatic_input(datasets_tab=datasets_tab, sample=wildcards.sample, library=wildcards.library)[1]
        # R1 = "data/reads/{sample}_{library}.R1.fastq.gz",
        # R2 = "data/reads/{sample}_{library}.R2.fastq.gz"
    output:
        out1P = "data/reads_filtered/{sample}_{library}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{library}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{library}_qc_U.fastq.gz",
    threads:
        config['read_processing']['trimmomatic']['threads']
    # version:
    #     subprocess.check_output("trimmomatic -version", shell=True)
    message:
        "Filtering read dataset {wildcards.sample}_{wildcards.library} with Trimmomatic. {wildcards}" # v{version}"
    log:
        log_dir + "/trimmomatic/{sample}_{library}_trimmomatic.log"
    #conda: "envs/environment.yaml"
    run:
        #trimmomatic_adapters_path = get_trimmomatic_adapters_path()
        shell("export tap=$(which trimmomatic | sed 's/bin\/trimmomatic/share\/trimmomatic\/adapters\/TruSeq3-PE.fa/g'); trimmomatic PE {params.options} -threads {threads} {input.R1} {input.R2} {params.out1P} {params.out1U} {params.out2P} {params.out2U} ILLUMINACLIP:$tap:2:30:10 {params.processing_options} &> {log}")
        shell("zcat {params.out1U} {params.out2U} | gzip > {output.out1U} && rm {params.out1U} {params.out2U}")

seq_type = "both"

rule map_MT_PE_SE:
    input:
        R1 = "data/reads_filtered/{sample}_{library}_qc_R1.fastq.gz",
        R2 = "data/reads_filtered/{sample}_{library}_qc_R2.fastq.gz",
        U = "data/reads_filtered/{sample}_{library}_qc_U.fastq.gz",
        gmap_db = genome_db(gmap_db_dir + "/{ref_organism}/{ref_organism}.done"),
        # gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    output:
        outmt_sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt.sam.gz"
    params:
        gmap_db_dir =  lambda wildcards: os.path.join(gmap_db_dir, wildcards.ref_organism),
        gmap_db = lambda wildcards: ref_organism_dict[wildcards.ref_organism].ref_genome_mt,
        RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample',
        uncompressed_output = lambda wildcards, output: output.outmt_sam.replace("_outmt.sam.gz", "_outmt.sam")
    log:
        log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/{sample}_{library}_{ref_organism}_map_MT_PE_SE.log"
    #conda: "envs/environment.yaml"
    threads:
        config["map"]["gmap_threads"]
    message: "Mapping reads for read dataset {wildcards.sample}_{wildcards.library} to {wildcards.ref_organism} mt genome"
    run:
        if seq_type == "pe":
            print("PE mode")
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} &> {log} && gzip {params.uncompressed_output} &>> {log}")
        if seq_type == "se":
            print("SE mode")
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} &> {log} && gzip {params.uncompressed_output} &>> {log}")
        elif seq_type == "both":
            print("PE + SE mode")
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} {input[2]} &> {log} && gzip {params.uncompressed_output} &>> {log}")

rule sam_to_ids:
    input:
        outmt_sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt.sam.gz"
    output:
        outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt1.ids",
        #outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt2.ids",
        outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt.ids",
        outmt_U1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt_U1.ids",
        outmt_U2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt_U2.ids",
    message:
        "Getting ids of mapped reads from {input.outmt_sam}"
    run:
        sam_to_ids(samfile=input.outmt_sam, outmt_PE=output.outmt1, outmt_U1=output.outmt_U1, outmt_U2=output.outmt_U2, 
                             outmt_SE=output.outmt, keep_orphans=True, return_dict=False, return_files=True)

rule ids_to_fastq_PE:
    input:
        outmt1 = rules.sam_to_ids.output.outmt1,# if config["keep_orphans"] \
                    #else rules.sam_to_ids.output.outmt1,
        R1 = rules.trimmomatic.output.out1P,
        R2 = rules.trimmomatic.output.out2P,
        #outmt = rules.sam_to_ids.output.outmt,
    output:
        outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt1.fastq.gz",
        outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt2.fastq.gz",
        #outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.fastq.gz",
    message:
        "Fetching reads in {input.outmt1}"
    params:
        out1_temp = lambda wildcards, output: output.outmt1.replace(".gz", ""),
        out2_temp = lambda wildcards, output: output.outmt2.replace(".gz", "") 
    run:
        run_seqtk_subset(seqfile=input.R1, id_list=input.outmt1, outseqfile=params.out1_temp)
        run_seqtk_subset(seqfile=input.R2, id_list=input.outmt1, outseqfile=params.out2_temp)
        shell("gzip {params.out1_temp}")
        shell("gzip {params.out2_temp}")

rule ids_to_fastq_SE:
    input:
        outmt = rules.sam_to_ids.output.outmt,# if config["keep_orphans"] \
                    #else rules.sam_to_ids.output.outmt,
        U1 = rules.trimmomatic.output.out1U,
        #U2 = rules.trimmomatic.output.out2U,
        #outmt = rules.sam_to_ids.output.outmt,
    output:
        outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt.fastq.gz",
        #outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt_PE_2.fastq.gz",
        #outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.fastq.gz",
    message:
        "Fetching reads in {input.outmt}"
    params:
        out1_temp = lambda wildcards, output: output.outmt1.replace(".gz", ""),
#        out2_temp = lambda wildcards, output: output.outmt2.replace(".gz", "") 
    run:
        run_seqtk_subset(seqfile=input.U1, id_list=input.outmt, outseqfile=params.out1_temp)
#        run_seqtk_subset(seqfile=input.R2, id_list=input.outmt1, outseqfile=params.out2_temp)
        shell("gzip {params.out1_temp}")
#        shell("gzip {params.out2_temp}")

rule ids_to_fastq_orphans:
    input:
        outmt_U1 = rules.sam_to_ids.output.outmt_U1,
        outmt_U2 = rules.sam_to_ids.output.outmt_U2,
        R1 = rules.trimmomatic.output.out1P,
        R2 = rules.trimmomatic.output.out2P,
        #outmt = rules.sam_to_ids.output.outmt,
    output:
        outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt_U1.fastq.gz",
        outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outmt_U2.fastq.gz",
        #outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.fastq.gz",
    message:
        "Fetching reads in {input.outmt_U1} and {input.outmt_U2}"
    params:
        out1_temp = lambda wildcards, output: output.outmt1.replace(".gz", ""),
        out2_temp = lambda wildcards, output: output.outmt2.replace(".gz", "") 
    run:
        run_seqtk_subset(seqfile=input.R1, id_list=input.outmt_U1, outseqfile=params.out1_temp)
        run_seqtk_subset(seqfile=input.R2, id_list=input.outmt_U2, outseqfile=params.out2_temp)
        shell("gzip {params.out1_temp}")
        shell("gzip {params.out2_temp}")

# rule map_nuclear_MT_SE:
#     input:
#         lambda wildcards: get_inputs_for_rule_map_nuclear_MT_SE(sample=wildcards.sample, library=wildcards.library,
#                                                                 ref_organism=wildcards.ref_organism,
#                                                                 keep_orphans=keep_orphans),
#         gmap_db = genome_db(gmap_db_dir + "/{ref_organism}/{ref_organism}.done")
#         #outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.fastq.gz",
#     output:
#         outS = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outS.sam.gz"
#     params:
#         gmap_db_dir =  lambda wildcards: os.path.join(gmap_db_dir, wildcards.ref_organism),
#         gmap_db =      lambda wildcards: "{ref_genome_mt}_{ref_genome_n}".format(ref_organism_dict[wildcards.ref_organism].ref_genome_mt,
#                                                                                  ref_organism_dict[wildcards.ref_organism].ref_genome_n),
#         uncompressed_output = lambda wildcards, output: output.outS.replace("_outS.sam.gz", "_outS.sam")
#     threads:
#         config["map"]["gmap_remap_threads"]
#     # version:
#     #     subprocess.getoutput(
#     #       gsnap --version
#     #       )
#     #     """
#     #conda: "envs/environment.yaml"
#     log:
#         logS = log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/{sample}_{library}_{ref_organism}_map_nuclear_MT_SE.log"
#     message:
#         "Mapping onto complete human genome (nuclear + mt)... SE reads"
#     run:
#         if os.path.isfile(input.outmt):
#             shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input[:-1]} &> {log.logS} && gzip {params.uncompressed_output} &>> {log.logS}")
#         else:
#             open(output.outS, 'a').close()
# 
# rule map_nuclear_MT_PE:
#     input:
#         gmap_db = genome_db(gmap_db_dir + "/{ref_organism}/{ref_organism}.done"),
#         outmt1 = rules.ids_to_fastq_PE.output.outmt1,
#         outmt2 = rules.ids_to_fastq_PE.output.outmt2,
#         # outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt_PE_1.fastq.gz",
#         # outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt_PE_2.fastq.gz",
#     output:
#         outP = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_outP.sam.gz"
#     params:
#         gmap_db_dir = config["map"]["gmap_db_dir"],
#         #gsnap_db_folder = config['map']['gsnap_db_folder'],
#         gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].replace(".chromosome", ""),
#         uncompressed_output = lambda wildcards, output: output.outP.replace("_outP.sam.gz", "_outP.sam")
#         #gsnap_db = config['map']['gsnap_n_mt_db']
#     threads:
#         config["map"]["gmap_remap_threads"]
#     # version:
#     #     subprocess.getoutput(
#     #       gsnap --version
#     #       )
#     #     """
#     #conda: "envs/environment.yaml"
#     log:
#         logP = log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/{sample}_{library}_{ref_organism}_map_nuclear_MT_PE.log"
#     message:
#         "Mapping onto complete human genome (nuclear + mt)... PE reads"
#     run:
#         if os.path.isfile(input.outmt1):
#             shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt1} {input.outmt2} &> {log.logP} && gzip {params.uncompressed_output} &>> {log.logP}")
#         else:
#             open(output.outP, 'a').close()

rule map_nuclear_MT_PE_SE:
    input:
        outmt1 = rules.ids_to_fastq_PE.output.outmt1,
        outmt2 = rules.ids_to_fastq_PE.output.outmt2,
        outmt_SE = lambda wildcards: get_inputs_for_rule_map_nuclear_MT_SE(sample=wildcards.sample, library=wildcards.library,
                                                                ref_organism=wildcards.ref_organism, keep_orphans=keep_orphans),
        gmap_db = genome_db(gmap_db_dir + "/{ref_organism}/{ref_organism}.done"),
    output:
        concordant_uniq = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_out_mt_n.concordant_uniq",
        concordant_circular = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_out_mt_n.concordant_circular",
        unpaired_uniq = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_out_mt_n.unpaired_uniq",
        unpaired_circular = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_out_mt_n.unpaired_circular",
        # outmt_n = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_out_mt_n.sam.gz"
    threads: 4
    log:
        log = log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/{sample}_{library}_{ref_organism}_map_nuclear_MT_PE_SE.log"
    params:
        gmap_db_dir =  lambda wildcards: os.path.join(gmap_db_dir, wildcards.ref_organism),
        gmap_db =      lambda wildcards: "{ref_genome_mt}_{ref_genome_n}".format(ref_genome_mt=ref_organism_dict[wildcards.ref_organism].ref_genome_mt,
                                                                                 ref_genome_n=ref_organism_dict[wildcards.ref_organism].ref_genome_n),
        # gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].replace(".chromosome", ""),
        out_basename = lambda wildcards, output: output.concordant_uniq.replace(".concordant_uniq", ""),
        RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample',
    run:
        input_files_PE = input[:2]   # PE files are always the first two
        input_files_SE = input[2:-1] # SE files are those from the third one to the one before the last one
        # create a single SE file to be compliant with gsnap
        input_files_SE_all = input_files_SE[0].replace("_outmt.fastq.gz", "_outmt_all.fastq.gz")
        if len(input_files_SE) > 1:
            with open(input_files_SE_all, 'wb') as input_files_SE_cat:
                for f in input_files_SE:
                    with open(f, 'rb') as rfp:
                        shutil.copyfileobj(rfp, input_files_SE_cat)
            input_files_SE = [input_files_SE_all]
        shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} --split-output={params.out_basename} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input_files_PE} {input_files_SE} &> {log}")#" && gzip {params.uncompressed_output} &>> {log}")
        # delete the single SE file
        os.remove(input_files_SE_all)

rule filtering_mt_alignments:
    input:
        rules.map_nuclear_MT_PE_SE.output,
#        mt_n_fasta = genome_db(rootdir + "/data/genomes/{ref_organism}_mt_n.fasta")
    output:
        sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT.sam.gz"
    message: "Filtering alignments in files {input}"
    run:
        ref_mt_fasta_header = lambda wildcards: get_mt_AC(genome_db_data=genome_db_data,
                                                            ref_organism=wildcards.ref_organism, reference_tab=reference_tab)
        filtering_report = cat_alignments(input, outfile=output.sam, ref_mt_fasta_header=ref_mt_fasta_header)
        print(filtering_report)

rule sam2bam:
    input:
        sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT.sam.gz",
    output:
        "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT.bam",
    message: "Converting {input.sam} to {output}"
    log: log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/sam2bam.log"
    #group: "variant_calling"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        zcat {input.sam} | samtools view -b -o {output} - &> {log}
        """

rule sort_bam:
    input:
        bam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT.bam"
    output:
        sorted_bam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT-sorted.bam"
    message: "Sorting {input.bam} to {output.sorted_bam}"
    params:
        TMP = check_tmp_dir(config["tmp_dir"])
    log: log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/sort_bam.log"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
    shell:
        """
        samtools sort -o {output.sorted_bam} -T {params.TMP} {input.bam} &> {log}
        """

rule mark_duplicates:
    input:
        sorted_bam = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT-sorted.bam"
    output:
        sorted_bam_md = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT-sorted.final.bam",
        metrics_file = "results/{sample}/map/OUT_{sample}_{library}_{ref_organism}/{sample}_{library}_{ref_organism}_OUT-sorted.final.metrics.txt"
    params:
        TMP = check_tmp_dir(config["tmp_dir"]),
        mark_duplicates = config["mark_duplicates"]
    message: "Removing duplicate reads from {input.sorted_bam}: {params.mark_duplicates}. Output: {output.sorted_bam_md}"
    log: log_dir + "/{sample}/OUT_{sample}_{library}_{ref_organism}/map/mark_duplicates.log"
    run:
        if params.mark_duplicates == True:
            shell("picard MarkDuplicates \
                INPUT={input.sorted_bam} \
                OUTPUT={output.sorted_bam_md} \
                METRICS_FILE={output.metrics_file} \
                ASSUME_SORTED=true \
                REMOVE_DUPLICATES=true \
                TMP_DIR={params.TMP}")
        else:
            shutil.copy2(input.sorted_bam, output.sorted_bam_md)
            with open(output.metrics_file, "w") as f:
                f.write("")

rule merge_bam:
    input:
        sorted_bams = lambda wildcards: get_sample_bamfiles(datasets_tab, res_dir="results",
                                                            sample=wildcards.sample,
                                                            ref_organism=wildcards.ref_organism)
                                                            # ref_genome_mt=wildcards.ref_genome_mt,
                                                            # ref_genome_n=wildcards.ref_genome_n)
    output:
        merged_bam = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.bam",
        merged_bam_index = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.bam.bai"
    log: log_dir + "/{sample}/{sample}_{ref_organism}_merge_bam.log"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        samtools merge {output.merged_bam} {input} &> {log}
        samtools index {output.merged_bam} {output.merged_bam_index}
        """

rule left_align_merged_bam:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.bam",
        merged_bam_index = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.bam.bai",
        mt_n_fasta = genome_db(rootdir + "/data/genomes/{ref_organism}_mt_n.fasta"),
        genome_index = genome_db(rootdir + "/data/genomes/{ref_organism}_mt_n.fasta.fai"),
        genome_dict = genome_db(rootdir + "/data/genomes/{ref_organism}.dict")
    output:
        merged_bam_left_realigned = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.bam"
    log: log_dir + "/{sample}/{sample}_{ref_organism}_left_align_merged_bam.log"
    params:
        source_dir = source_dir
    message: "Realigning indels in {input.merged_bam} with GATK 3.6 - LeftAlignIndels"
    shell:
        """
        gatk-framework -Xmx6G \
            -R {input.mt_n_fasta} \
            -T LeftAlignIndels \
            -I {input.merged_bam} \
            -o {output.merged_bam_left_realigned} \
            --filter_reads_with_N_cigar
        """

rule clip_bam:
    input:
        merged_bam_left_realigned = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.bam"
    output:
        merged_bam_left_realigned_clipped = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.clip.bam"
    log: log_dir + "/{sample}/{sample}_{ref_organism}_left_align_merged_bam_clip.log"
    message: "Clipping alignment ends {input.merged_bam_left_realigned} with bam trimBam"
    run:
        shell("bam trimBam \
            {input.merged_bam_left_realigned} \
            {output.merged_bam_left_realigned_clipped} \
            10 \
            -c")

rule clip_bam_recalculate_MD:
    input:
        merged_bam_left_realigned_clipped = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.clip.bam"
    output:
        merged_bam_left_realigned_clipped_newMD = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.clip.calmd.bam"
    params:
        ref_mt_fasta = lambda wildcards: rootdir + "/data/genomes/{ref_genome_mt_file}".format(
            ref_genome_mt_file=get_mt_fasta(reference_tab, ref_organism=wildcards.ref_organism, field="ref_genome_mt_file")
            ),
        TMP = check_tmp_dir(config["tmp_dir"])
    log: log_dir + "/{sample}/{sample}_{ref_organism}_left_align_merged_bam_clip_calmd.log"
    message: "Recalculating MD string for clipped alignments {input.merged_bam_left_realigned_clipped} with bam trimBam"
    run:
        shell("samtools calmd -b \
            {input.merged_bam_left_realigned_clipped} \
            {params.ref_mt_fasta} 2> {log} 1> {params.TMP}/$(basename {output.merged_bam_left_realigned_clipped_newMD})")
        shell("samtools sort -o {output.merged_bam_left_realigned_clipped_newMD} \
            -T {params.TMP} \
            {params.TMP}/$(basename {output.merged_bam_left_realigned_clipped_newMD}) &>> {log}")

rule bam2cov:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.clip.calmd.bam" if config["trimBam"] \
                        else "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.bam",
        #genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
    output:
        bam_cov = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.bam.cov"
    # params:
    #     genome_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    params:
        quality = config['mtvcf_main_analysis']['Q'],
    message: "Generating coverage file {output.bam_cov} from {input.merged_bam}"
    log: log_dir + "/{sample}/{sample}_{ref_organism}_bam2cov.log"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
    shell:
        """
        samtools depth \
        -d 0 \
        -q {params.quality} \
        -a {input.merged_bam} > {output.bam_cov}
        """

rule make_single_VCF:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.clip.calmd.bam" if config["trimBam"] \
                        else "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.bam",
        bam_cov = "results/{sample}/map/{sample}_{ref_organism}_OUT-sorted.realign.bam.cov",
    output:
        single_vcf = "results/{sample}/{sample}_{ref_organism}.vcf.gz",
        single_bed = "results/{sample}/{sample}_{ref_organism}.bed",
        single_fasta = "results/{sample}/{sample}_{ref_organism}.fasta"
    params:
        ref_mt_fasta = lambda wildcards: rootdir + "/data/genomes/{ref_genome_mt_file}".format(
            ref_genome_mt_file=get_mt_fasta(reference_tab, ref_organism=wildcards.ref_organism, field="ref_genome_mt_file")
            ),
        TMP = check_tmp_dir(config["tmp_dir"]),
        tail = config['mtvcf_main_analysis']['tail'],
        quality = config['mtvcf_main_analysis']['Q'],
        minrd = config['mtvcf_main_analysis']['minrd'],
        tail_mismatch = config['mtvcf_main_analysis']['tail_mismatch']
    message: "Processing {input.merged_bam} to get VCF {output.single_vcf}"
    run:
        # function (and related ones) from mtVariantCaller
        tmp_sam = os.path.split(input.merged_bam)[1].replace(".bam", ".sam")
        shell("samtools view {merged_bam} > {tmp_dir}/{tmp_sam}".format(merged_bam=input.merged_bam,
                                                                        tmp_dir=params.TMP,
                                                                        tmp_sam=tmp_sam))

        vcf_dict = mtvcf_main_analysis(sam_file="{tmp_dir}/{tmp_sam}".format(tmp_dir=params.TMP, tmp_sam=tmp_sam),
                                       coverage_data_file=input.bam_cov, name2=wildcards.sample,
                                       tail=params.tail, Q=params.quality, minrd=params.minrd,
                                       ref_mt=params.ref_mt_fasta, tail_mismatch=params.tail_mismatch)
        # ref_genome_mt will be used in the VCF descriptive field
        # seq_name in the VCF data
        seq_name, seq_length = get_seq_name(params.ref_mt_fasta)
        VCF_RECORDS = VCFoutput(vcf_dict, reference=seq_name, seq_name=seq_name, seq_length=seq_length,
                                vcffile=output.single_vcf)
        bed_output(VCF_RECORDS, seq_name=seq_name, bedfile=output.single_bed)
        # fasta output
        #contigs = pileup2mt_table(pileup=input.pileup, fasta=params.ref_mt_fasta, mt_table=in.mt_table)
        #mt_table_data = pileup2mt_table(pileup=input.pileup, ref_fasta=params.ref_mt_fasta)
        gapped_fasta = sam_cov_handle2gapped_fasta(coverage_data_file=input.bam_cov,
                                                   ref_mt=params.ref_mt_fasta)
        contigs = gapped_fasta2contigs(gapped_fasta=gapped_fasta)
        fasta_output(vcf_dict=vcf_dict, ref_mt=params.ref_mt_fasta, fasta_out=output.single_fasta,
                    contigs=contigs)

rule index_VCF:
    input:
        single_vcf = "results/{sample}/{sample}_{ref_organism}.vcf.gz",
    output:
        index_vcf = "results/{sample}/{sample}_{ref_organism}.vcf.gz.csi"
    #conda: "envs/bcftools.yaml"
    message: "Compressing and indexing {input.single_vcf}"
    run:
        shell("bcftools index {input.single_vcf}")

rule merge_VCF:
    input:
        single_vcf_list = lambda wildcards: get_genome_single_vcf_files(analysis_tab,
                                                                        ref_organism=wildcards.ref_organism),
        index_vcf = lambda wildcards: get_genome_single_vcf_index_files(analysis_tab,
                                                                        ref_organism=wildcards.ref_organism),
        #single_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz",
        #index_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz.csi"
    output:
        merged_vcf = "results/vcf/{ref_organism}.vcf"
    message: "Merging vcf files for reference organism: {wildcards.ref_organism}"
    #conda: "envs/bcftools.yaml"
    run:
        if len(input.single_vcf_list) == 1:
            shutil.copy2(input.single_vcf_list[0], output.merged_vcf+".gz")
            shell("gunzip {output.merged_vcf}.gz")
        else:
            shell("bcftools merge {input.single_vcf_list} -O v -o {output.merged_vcf}")
