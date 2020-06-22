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

from Bio import SeqIO
import numpy as np
import pandas as pd
from sqlalchemy import create_engine

for path in sys.path:
    if "snakefiles" in path:
        sys.path.append(path.replace("/snakefiles", ""))

from modules.BEDoutput import bed_output, fasta_output
from modules.config_parsers import (
    fastqc_outputs, get_bed_files, get_datasets_for_symlinks,
    get_fasta_files, get_genome_files, get_genome_single_vcf_files,
    get_genome_single_vcf_index_files, get_genome_vcf_files, get_mt_genomes, get_mt_fasta,
    get_sample_bamfiles, get_symlinks
)
from modules.filter_alignments import filter_alignments
from modules.general import (
    check_tmp_dir, gapped_fasta2contigs, get_seq_name, sam_to_fastq, sam_cov_handle2gapped_fasta
)
from modules.mtVariantCaller import mtvcf_main_analysis, VCFoutput

source_dir = Path(os.path.dirname(workflow.snakefile)).parent
#source_dir = os.path.abspath(os.path.join(".", os.pardir))
#localrules: bam2pileup, index_genome, pileup2mt_table, make_single_VCF
localrules: index_genome, merge_VCF, index_VCF, dict_genome, symlink_libraries

# fields: sample  ref_genome_mt   ref_genome_n
analysis_tab = pd.read_table("data/analysis.tab", sep = "\t", comment='#')
reference_tab = (pd.read_table("data/reference_genomes.tab", sep = "\t", comment='#')
                 .set_index("ref_genome_mt", drop=False))
datasets_tab = pd.read_table("data/datasets.tab", sep = "\t", comment='#')

configfile: "config.yaml"
res_dir = config["results"]
map_dir = config["map_dir"]
log_dir = config["log_dir"]
gmap_db_dir = config["map"]["gmap_db_dir"]

wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in list(set(analysis_tab['sample']))]),
    ref_genome_mt = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_mt']))]),
    ref_genome_n = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_n']))])

outpaths = get_mt_genomes(analysis_tab)

target_inputs = [
    outpaths ]

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

rule fastqc_raw:
    input:
        R1 = "data/reads/{sample}_{library}.R1.fastq.gz",
        R2 = "data/reads/{sample}_{library}.R2.fastq.gz",
        # R1 = "data/reads/{dataset_basename}_R1_001.fastq.gz",
        # R2 = "data/reads/{dataset_basename}_R2_001.fastq.gz"
    output:
        html_report_R1 = "results/fastqc_raw/{sample}_{library}.R1_fastqc.html",
        html_report_R2 = "results/fastqc_raw/{sample}_{library}.R2_fastqc.html",
        # html_report_R1 = "results/fastqc_raw/{dataset_basename}_R1_001_fastqc.html",
        # html_report_R2 = "results/fastqc_raw/{dataset_basename}_R2_001_fastqc.html",
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

rule make_mt_gmap_db:
    input:
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}",
                                                   ref_genome_mt_file=get_genome_files(reference_tab,
                                                                                       wildcards.ref_genome_mt,
                                                                                       "ref_genome_mt_file"))
    output:
        gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].replace(".chromosome", "")
    message: "Generating gmap db for mt genome: {input.mt_genome_fasta}.\nWildcards: {wildcards}"
    log: "logs/gmap_build/{ref_genome_mt}.log"
    #conda: "envs/environment.yaml"
    shell:
        """
        #module load gsnap
        gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -s none {input.mt_genome_fasta} &> {log}
        """

rule make_mt_n_gmap_db:
    input:
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}",
                                                   ref_genome_mt_file=get_genome_files(reference_tab,
                                                                                       wildcards.ref_genome_mt,
                                                                                       "ref_genome_mt_file")),
        n_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_n_file}",
                                                  ref_genome_n_file=get_genome_files(reference_tab,
                                                                                     wildcards.ref_genome_mt,
                                                                                     "ref_genome_n_file"))
    output:
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome",
        mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        # gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].split(".")[0]
        gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].replace(".chromosome", "")
    message: "Generating gmap db for mt + n genome: {input.mt_genome_fasta},{input.n_genome_fasta}"
    log: "logs/gmap_build/{ref_genome_mt}_{ref_genome_n}.log"
    #conda: "envs/environment.yaml"
    shell:
        """
        cat {input.mt_genome_fasta} {input.n_genome_fasta} > {output.mt_n_fasta}
        gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -s none {output.mt_n_fasta} &> {log}
        # rm {input.mt_genome_fasta}_{input.n_genome_fasta}.fasta
        """

rule fastqc_filtered:
    input:
        out1P = "data/reads_filtered/{sample}_{library}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{library}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{library}_qc_U.fastq.gz",
    output:
        html_report_R1 = "results/fastqc_filtered/{sample}_{library}.R1_fastqc.html",
        html_report_R2 = "results/fastqc_filtered/{sample}_{library}.R2_fastqc.html",
        html_report_U = "results/fastqc_filtered/{sample}_{library}.U_fastqc.html",
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
        R1 = "data/reads/{sample}_{library}.R1.fastq.gz",
        R2 = "data/reads/{sample}_{library}.R2.fastq.gz"
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
        gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    output:
        outmt_sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        gmap_db = lambda wildcards: wildcards.ref_genome_mt,
        RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample',
        uncompressed_output = lambda wildcards, output: output.outmt_sam.replace("_outmt.sam.gz", "_outmt.sam")
    log:
        log_dir + "/{sample}/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{library}_{ref_genome_mt}_map_MT_PE_SE.log"
    #conda: "envs/environment.yaml"
    threads:
        config["map"]["gmap_threads"]
    message: "Mapping reads for read dataset {wildcards.sample}_{wildcards.library} to {wildcards.ref_genome_mt} mt genome"
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

rule sam2fastq:
    input:
        outmt_sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.sam.gz"
        #outmt_sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_outmt.sam.gz"
    output:
        outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt1.fastq.gz",
        outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt2.fastq.gz",
        outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.fastq.gz",
        #log = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/sam2fastq.done"
    #conda: "envs/environment.yaml"
    message:
        "Converting {input.outmt_sam} to FASTQ"
    run:
        sclipped = sam_to_fastq(samfile=input.outmt_sam, outmt1=output.outmt1,
                             outmt2=output.outmt2, outmt=output.outmt, do_softclipping=True)
        print("{} reads with soft-clipping > 1/3 of their length".format(sclipped))

rule map_nuclear_MT_SE:
    input:
        outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.fastq.gz",
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome"
    output:
        outS = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_outS.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].replace(".chromosome", ""),
        uncompressed_output = lambda wildcards, output: output.outS.replace("_outS.sam.gz", "_outS.sam")
    threads:
        config["map"]["gmap_remap_threads"]
    # version:
    #     subprocess.getoutput(
    #       gsnap --version
    #       )
    #     """
    #conda: "envs/environment.yaml"
    log:
        logS = log_dir + "/{sample}/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_map_nuclear_MT_SE.log"
    message:
        "Mapping onto complete human genome (nuclear + mt)... SE reads"
    run:
        if os.path.isfile(input.outmt):
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt} &> {log.logS} && gzip {params.uncompressed_output} &>> {log.logS}")
        else:
            open(output.outS, 'a').close()

rule map_nuclear_MT_PE:
    input:
        outmt1 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt1.fastq.gz",
        outmt2 = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt2.fastq.gz",
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome"
    output:
        outP = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_outP.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].replace(".chromosome", ""),
        uncompressed_output = lambda wildcards, output: output.outP.replace("_outP.sam.gz", "_outP.sam")
        #gsnap_db = config['map']['gsnap_n_mt_db']
    threads:
        config["map"]["gmap_remap_threads"]
    # version:
    #     subprocess.getoutput(
    #       gsnap --version
    #       )
    #     """
    #conda: "envs/environment.yaml"
    log:
        logP = log_dir + "/{sample}/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_map_nuclear_MT_PE.log"
    message:
        "Mapping onto complete human genome (nuclear + mt)... PE reads"
    run:
        if os.path.isfile(input.outmt1):
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt1} {input.outmt2} &> {log.logP} && gzip {params.uncompressed_output} &>> {log.logP}")
        else:
            open(output.outP, 'a').close()

rule filtering_mt_alignments:
    input:
        outmt = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_outmt.sam.gz",
        outS = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_outS.sam.gz",
        outP = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_outP.sam.gz"
    output:
        sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(
            ref_genome_mt_file=get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file")
        )
    #conda: "envs/environment.yaml"
    threads: 1
    message: "Filtering alignments in file {input.outmt} by checking alignments in {input.outS} and {input.outP}"
    run:
        filter_alignments(outmt=input.outmt, outS=input.outS, outP=input.outP, OUT=output.sam,
                          ref_mt_fasta=params.ref_mt_fasta)

rule sam2bam:
    input:
        sam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz",
    output:
        "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT.bam",
    message: "Converting {input.sam} to {output}"
    log: log_dir + "/{sample}/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/map/sam2bam.log"
    #group: "variant_calling"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        zcat {input.sam} | samtools view -b -o {output} - &> {log}
        """

rule sort_bam:
    input:
        bam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT.bam"
    output:
        sorted_bam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam"
    message: "Sorting {input.bam} to {output.sorted_bam}"
    params:
        TMP = check_tmp_dir(config["tmp_dir"])
    log: log_dir + "/{sample}/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/map/sort_bam.log"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
    shell:
        """
        samtools sort -o {output.sorted_bam} -T {params.TMP} {input.bam} &> {log}
        # samtools sort -o {output.sorted_bam} -T ${{TMP}} {input.bam}
        """

rule mark_duplicates:
    input:
        sorted_bam = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam"
    output:
        sorted_bam_md = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.final.bam",
        metrics_file = "results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.final.metrics.txt"
    params:
        TMP = check_tmp_dir(config["tmp_dir"]),
        mark_duplicates = config["mark_duplicates"]
    message: "Removing duplicate reads from {input.sorted_bam}: {params.mark_duplicates}. Output: {output.sorted_bam_md}"
    log: log_dir + "/{sample}/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/map/mark_duplicates.log"
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
                                                            ref_genome_mt=wildcards.ref_genome_mt,
                                                            ref_genome_n=wildcards.ref_genome_n)
    output:
        merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam",
        merged_bam_index = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam.bai"
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_merge_bam.log"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        samtools merge {output.merged_bam} {input} &> {log}
        samtools index {output.merged_bam} {output.merged_bam_index}
        """

rule index_genome:
    input:
        mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
    message: "Indexing {input.mt_n_fasta} with samtools faidx"
    log: log_dir + "/{ref_genome_mt}_{ref_genome_n}.samtools_index.log"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        samtools faidx {input.mt_n_fasta} &> {log}
        """

rule dict_genome:
    input:
        mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        genome_dict = "data/genomes/{ref_genome_mt}_{ref_genome_n}.dict"
    message: "Creating .dict of {input.mt_n_fasta} with picard CreateSequenceDictionary"
    log: log_dir + "/{ref_genome_mt}_{ref_genome_n}.picard_dict.log"
    #conda: "envs/samtools_biopython.yaml"
    run:
        shell("picard CreateSequenceDictionary R={input.mt_n_fasta} O={output.genome_dict}")

rule left_align_merged_bam:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam",
        merged_bam_index = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam.bai",
        mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta",
        genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai",
        genome_dict = "data/genomes/{ref_genome_mt}_{ref_genome_n}.dict"
    output:
        merged_bam_left_realigned = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam"
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_left_align_merged_bam.log"
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
        merged_bam_left_realigned = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam"
    output:
        merged_bam_left_realigned_clipped = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.clip.bam"
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_left_align_merged_bam_clip.log"
    message: "Clipping alignment ends {input.merged_bam_left_realigned} with bam trimBam"
    run:
        shell("bam trimBam \
            {input.merged_bam_left_realigned} \
            {output.merged_bam_left_realigned_clipped} \
            10 \
            -c")

rule clip_bam_recalculate_MD:
    input:
        merged_bam_left_realigned_clipped = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.clip.bam"
    output:
        merged_bam_left_realigned_clipped_newMD = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.clip.calmd.bam"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(
            ref_genome_mt_file=get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file")
        ),
        TMP = check_tmp_dir(config["tmp_dir"])
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_left_align_merged_bam_clip_calmd.log"
    message: "Recalculating MD string for clipped alignments {input.merged_bam_left_realigned_clipped} with bam trimBam"
    run:
        shell("samtools calmd -b \
            {input.merged_bam_left_realigned_clipped} \
            {params.ref_mt_fasta} 2> {log} 1> {params.TMP}/$(basename {output.merged_bam_left_realigned_clipped_newMD})")
        shell("samtools sort -o {output.merged_bam_left_realigned_clipped_newMD} \
            -T {params.TMP} \
            {params.TMP}/$(basename {output.merged_bam_left_realigned_clipped_newMD}) &>> {log}")

# rule bam2pileup:
#     input:
#         merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam",
#         genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
#     output:
#         pileup = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
#     params:
#         genome_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
#     message: "Generating pileup {output.pileup} from {input.merged_bam}"
#     log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_bam2pileup.log"
#     #conda: "envs/samtools_biopython.yaml"
#     #group: "variant_calling"
#     shell:
#         """
#         samtools mpileup -B -f {params.genome_fasta} -o {output.pileup} {input.merged_bam} &> {log}
#         """
#
# rule pileup2mt_table:
#     input:
#         pileup = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
#         # pileup = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
#     output:
#         mt_table = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-mt_table.txt"
#     params:
#         ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
#     message: "Generating mt_table {output.mt_table} from {input.pileup}, ref mt: {params.ref_mt_fasta}"
#     #conda: "envs/environment.yaml"
#     #group: "variant_calling"
#     run:
#         mt_table_data = pileup2mt_table(pileup=input.pileup, ref_fasta=params.ref_mt_fasta)
#         write_mt_table(mt_table_data=mt_table_data, mt_table_file=output.mt_table)

#         name="some/file.txt" if config["condition"] else "other/file.txt"
rule bam2cov:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.clip.calmd.bam" if config["trimBam"] \
                        else "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam",
        #genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
    output:
        bam_cov = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam.cov"
    # params:
    #     genome_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    params:
        quality = config['mtvcf_main_analysis']['Q'],
    message: "Generating coverage file {output.bam_cov} from {input.merged_bam}"
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_bam2cov.log"
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
        merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.clip.calmd.bam" if config["trimBam"] \
                        else "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam",
        bam_cov = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam.cov",
        # sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz",
        #mt_table = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-mt_table.txt",
        #pileup = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
        # sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz",
        # mt_table = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-mt_table.txt"
    output:
        single_vcf = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz",
        single_bed = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.bed",
        single_fasta = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.fasta"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(
            ref_genome_mt_file=get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file")
            ),
        TMP = check_tmp_dir(config["tmp_dir"]),
        tail = config['mtvcf_main_analysis']['tail'],
        quality = config['mtvcf_main_analysis']['Q'],
        minrd = config['mtvcf_main_analysis']['minrd'],
        tail_mismatch = config['mtvcf_main_analysis']['tail_mismatch']
    message: "Processing {input.merged_bam} to get VCF {output.single_vcf}"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
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
        seq_name = get_seq_name(params.ref_mt_fasta)
        VCF_RECORDS = VCFoutput(vcf_dict, reference=wildcards.ref_genome_mt, seq_name=seq_name,
                                vcffile=output.single_vcf)
        BEDoutput(VCF_RECORDS, seq_name=seq_name, bedfile=output.single_bed)
        # fasta output
        #contigs = pileup2mt_table(pileup=input.pileup, fasta=params.ref_mt_fasta, mt_table=in.mt_table)
        #mt_table_data = pileup2mt_table(pileup=input.pileup, ref_fasta=params.ref_mt_fasta)
        gapped_fasta = sam_cov_handle2gapped_fasta(sam_cov_data=sam_cov_dict,
                                                   ref_mt=params.ref_mt_fasta)
        contigs = gapped_fasta2contigs(gapped_fasta=gapped_fasta)
        FASTAoutput(vcf_dict=vcf_dict, ref_mt=params.ref_mt_fasta, fasta_out=output.single_fasta,
                    contigs=contigs)

rule index_VCF:
    input:
        single_vcf = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz",
    output:
        index_vcf = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz.csi"
    #conda: "envs/bcftools.yaml"
    message: "Compressing and indexing {input.single_vcf}"
    run:
        shell("bcftools index {input.single_vcf}")

rule merge_VCF:
    input:
        single_vcf_list = lambda wildcards: get_genome_single_vcf_files(analysis_tab,
                                                                        ref_genome_mt=wildcards.ref_genome_mt),
        index_vcf = lambda wildcards: get_genome_single_vcf_index_files(analysis_tab,
                                                                        ref_genome_mt=wildcards.ref_genome_mt),
        #single_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz",
        #index_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz.csi"
    output:
        merged_vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.vcf"
    message: "Merging vcf files for mt reference genome: {wildcards.ref_genome_mt}"
    #conda: "envs/bcftools.yaml"
    run:
        shell("bcftools merge {input.single_vcf_list} -O v -o {output.merged_vcf}")
