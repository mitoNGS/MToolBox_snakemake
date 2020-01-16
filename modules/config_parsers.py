#!/usr/bin/env python
import os
from typing import List

import pandas as pd
from snakemake.io import expand


# TODO: infolder is not used anywhere
def get_datasets_for_symlinks(df, sample=None, library=None, d=None,
                              infolder="data/reads", outfolder="data/reads"):
    dataset_file = None
    for row in df.itertuples():
        if (getattr(row, "sample") == sample and
                getattr(row, "library") == int(library)):
            dataset_file = os.path.join(outfolder, getattr(row, d))

    return dataset_file


# TODO: infolder is not used anywhere
def get_symlinks(df, analysis_tab=None,
                 infolder="data/reads", outfolder="data/reads"):
    outpaths = []
    # TODO: convert .iterrows() to .itertuples() for efficiency
    for i, l in df.iterrows():
        if l["sample"] in list(analysis_tab["sample"]):
            outpaths.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}.R1.fastq.gz".format(
                        sample=l["sample"], library=l["library"])
                )
            )
            outpaths.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}.R2.fastq.gz".format(
                        sample=l["sample"], library=l["library"])
                )
            )
    return outpaths


def get_genome_single_vcf_files(df, res_dir="results", ref_genome_mt=None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
            outpaths.append(("{results}/{sample}/{sample}_{ref_genome_mt}_"
                             "{ref_genome_n}.vcf.gz").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=getattr(row, "ref_genome_mt"),
                ref_genome_n=getattr(row, "ref_genome_n")))

    return list(set(outpaths))


# TODO: library is not used anywhere
def get_sample_bamfiles(df, res_dir="results", sample=None, library=None,
                        ref_genome_mt=None, ref_genome_n=None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "sample") == sample:
            bam_file = ("{sample}_{library}_{ref_genome_mt}_"
                        "{ref_genome_n}_OUT-sorted.final.bam").format(
                sample=sample,
                library=getattr(row, "library"),
                ref_genome_mt=ref_genome_mt,
                ref_genome_n=ref_genome_n)
            out_folder = "OUT_{base}".format(
                base=bam_file.replace("_OUT-sorted.final.bam", ""))
            outpaths.append("{results}/{sample}/map/{out_folder}/{bam_file}".format(
                results=res_dir,
                bam_file=bam_file,
                sample=sample,
                out_folder=out_folder))

    return outpaths


def get_genome_single_vcf_index_files(df, res_dir="results", ref_genome_mt=None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
            outpaths.append(
                ("{results}/{sample}/{sample}_{ref_genome_mt}_"
                 "{ref_genome_n}.vcf.gz.csi").format(
                    results=res_dir,
                    sample=getattr(row, "sample"),
                    ref_genome_mt=getattr(row, "ref_genome_mt"),
                    ref_genome_n=getattr(row, "ref_genome_n")))

    return list(set(outpaths))


def get_genome_vcf_files(df, res_dir="results/vcf"):
    outpaths = set()
    for row in df.itertuples():
        outpaths.add("{results}/{ref_genome_mt}_{ref_genome_n}.vcf".format(
            results=res_dir,
            ref_genome_mt=getattr(row, "ref_genome_mt"),
            ref_genome_n=getattr(row, "ref_genome_n")))
    outpaths = list(outpaths)
    return outpaths


def get_bed_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append(
            ("{results}/{sample}/{sample}_"
             "{ref_genome_mt}_{ref_genome_n}.bed").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=getattr(row, "ref_genome_mt"),
                ref_genome_n=getattr(row, "ref_genome_n")))

    return outpaths


def get_fasta_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append(
            ("{results}/{sample}/{sample}_"
             "{ref_genome_mt}_{ref_genome_n}.fasta").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=getattr(row, "ref_genome_mt"),
                ref_genome_n=getattr(row, "ref_genome_n")))

    return outpaths


def get_haplo_prediction_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append(
            ("{results}/{sample}/{sample}_"
             "{ref_genome_mt}_{ref_genome_n}.csv").format(
                results=res_dir,
                sample=getattr(row, "sample"),
                ref_genome_mt=getattr(row, "ref_genome_mt"),
                ref_genome_n=getattr(row, "ref_genome_n")))
    return outpaths


def get_genome_files(df, ref_genome_mt, field):
    return expand(df.loc[ref_genome_mt, field])


def get_mt_genomes(df: pd.DataFrame) -> List[str]:
    """ Return a list of unique mt genome identifiers from the
    given dataframe.

    Args:
        df: input pandas DataFrame

    Returns:
        list of str
    """
    return df["ref_genome_mt"].unique().tolist()


def get_mt_fasta(df, ref_genome_mt, field):
    return df.loc[df['ref_genome_mt'] == ref_genome_mt, field][0]


# TODO: infolder and ext are not used anywhere
def fastqc_raw_outputs(datasets_tab: pd.DataFrame,
                       analysis_tab: pd.DataFrame,
                       infolder="data/reads",
                       outfolder: str = "results/fastqc_raw",
                       ext=".fastq.gz") -> List[str]:
    """ Return a list of output filenames where FastQC results
    will be stored.

    Args:
        datasets_tab: input pandas DataFrame with fastq filenames
        analysis_tab: input pandas DataFrame with analysis details
        outfolder: output directory where FastQC results will be
            saved

    Returns:
        list of str paths
    """
    fastqc_out = []
    for _, l in datasets_tab.iterrows():
        if l["sample"] in analysis_tab["sample"].tolist():
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}.R1_fastqc.html".format(
                        sample=l["sample"], library=l["library"])
                )
            )
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}.R2_fastqc.html".format(
                        sample=l["sample"], library=l["library"])
                )
            )
    return fastqc_out


# TODO: infolder and ext are not used anywhere
def fastqc_filtered_outputs(datasets_tab, analysis_tab=None,
                            infolder="data/reads",
                            outfolder="results/fastqc_filtered",
                            ext="_001.fastq.gz"):
    fastqc_out = []
    for i, l in datasets_tab.iterrows():
        if l["sample"] in list(analysis_tab["sample"]):
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}_qc_R1_fastqc.html".format(
                        sample=l["sample"], library=l["library"])
                )
            )
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}_qc_R2_fastqc.html".format(
                        sample=l["sample"], library=l["library"])))
            fastqc_out.append(
                os.path.join(
                    outfolder,
                    "{sample}_{library}_qc_U_fastqc.html".format(
                        sample=l["sample"], library=l["library"])))
    return fastqc_out
