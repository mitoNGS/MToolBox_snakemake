#!/usr/bin/env python
import pandas as pd

def get_genome_single_vcf_files(df, res_dir="results", ref_genome_mt = None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
            outpaths.append("{results}/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz".format(results = res_dir, \
                                                                                                                        sample = getattr(row, "sample"), \
                                                                                                                        ref_genome_mt = getattr(row, "ref_genome_mt"), \
                                                                                                                        ref_genome_n = getattr(row, "ref_genome_n")))
            # outpaths.append("{results}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz".format(results = res_dir, \
            #                                                                                                             sample = getattr(row, "sample"), \
            #                                                                                                             ref_genome_mt = getattr(row, "ref_genome_mt"), \
            #                                                                                                             ref_genome_n = getattr(row, "ref_genome_n")))
    return outpaths

def get_sample_bamfiles(df, res_dir="results", sample = None, ref_genome_mt = None, ref_genome_n = None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "sample") == sample:
            bam_file = getattr(row, "R1").replace("_R1_001.fastq.gz", "_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.final.bam".format(ref_genome_mt = ref_genome_mt, ref_genome_n = ref_genome_n))
            out_folder = "OUT_{base}".format(base = bam_file.replace("_OUT-sorted.final.bam", ""))
            outpaths.append("{results}/{sample}/map/{out_folder}/{bam_file}".format(results = res_dir, bam_file = bam_file, sample = sample, out_folder = out_folder))
            # outpaths.append("{results}/{sample}/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam".format(results = res_dir, \
            #                                                                                                                                                                                     sample = ))
            #
            # outpaths.append("{results}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz".format(results = res_dir, \
            #                                                                                                             sample = getattr(row, "sample"), \
            #                                                                                                             ref_genome_mt = getattr(row, "ref_genome_mt"), \
            #                                                                                                             ref_genome_n = getattr(row, "ref_genome_n")))
    return outpaths

def get_genome_single_vcf_index_files(df, res_dir="results", ref_genome_mt = None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
            outpaths.append("{results}/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz.csi".format(results = res_dir, \
                                                                                                                        sample = getattr(row, "sample"), \
                                                                                                                        ref_genome_mt = getattr(row, "ref_genome_mt"), \
                                                                                                                        ref_genome_n = getattr(row, "ref_genome_n")))
            # outpaths.append("{results}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz.csi".format(results = res_dir, \
            #                                                                                                             sample = getattr(row, "sample"), \
            #                                                                                                             ref_genome_mt = getattr(row, "ref_genome_mt"), \
            #                                                                                                             ref_genome_n = getattr(row, "ref_genome_n")))
    return outpaths

def get_genome_vcf_files(df, res_dir="results/vcf"):
    outpaths = set()
    for row in df.itertuples():
        outpaths.add("{results}/{ref_genome_mt}_{ref_genome_n}.vcf".format(results = res_dir, \
                                                                            ref_genome_mt = getattr(row, "ref_genome_mt"), \
                                                                            ref_genome_n = getattr(row, "ref_genome_n")))
    outpaths = list(outpaths)
    return outpaths

def get_bed_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{results}/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.bed".format(results = res_dir, \
                                                                                                                    sample = getattr(row, "sample"), \
                                                                                                                    ref_genome_mt = getattr(row, "ref_genome_mt"), \
                                                                                                                    ref_genome_n = getattr(row, "ref_genome_n")))
        # outpaths.append("{results}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.bed".format(results = res_dir, \
        #                                                                                                             sample = getattr(row, "sample"), \
        #                                                                                                             ref_genome_mt = getattr(row, "ref_genome_mt"), \
        #                                                                                                             ref_genome_n = getattr(row, "ref_genome_n")))
    return outpaths

def get_genome_files(df, ref_genome_mt, field):
    return expand(df.loc[ref_genome_mt, field])

def get_mt_genomes(df):
    return list(set(df['ref_genome_mt']))

def get_mt_fasta(df, ref_genome_mt, field):
    return df.loc[df['ref_genome_mt'] == ref_genome_mt, field][0]