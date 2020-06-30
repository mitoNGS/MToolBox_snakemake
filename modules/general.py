#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import bz2
import gzip
import os
import re
import resource
import sys
from typing import Union
from types import SimpleNamespace
from snakemake import shell

from Bio import SeqIO
from Bio.Seq import reverse_complement

from modules.constants import CLEV, COV, DIUPAC, GLEN, MQUAL

shell.prefix("set -euo pipefail;")

def is_compr_file(f):
    with gzip.open(f, 'r') as fh:
        try:
            fh.read(1)
            return True
        except OSError:
            return False

def trimmomatic_input(datasets_tab=None, sample=None, library=None):
    """Input is supposed to be
    
    data/reads/{sample}_{library}.R1.fastq.gz"""
    outpaths = []
    for i, l in datasets_tab.iterrows():
        if l["sample"] == sample and str(l["library"]) == str(library):
            #print(l)
            if is_compr_file("data/reads/{}".format(l["R1"])):
                outpaths.append("data/reads/{sample}_{library}.R1.fastq.gz".format(sample=l["sample"], library=l["library"]))
                outpaths.append("data/reads/{sample}_{library}.R2.fastq.gz".format(sample=l["sample"], library=l["library"]))
            else:
                outpaths.append("data/reads/{sample}_{library}.R1.fastq".format(sample=l["sample"], library=l["library"]))
                outpaths.append("data/reads/{sample}_{library}.R2.fastq".format(sample=l["sample"], library=l["library"]))
    return outpaths

# def get_symlinks(df, analysis_tab=None,
#                  infolder="data/reads", outfolder="data/reads"):
#     outpaths = []
#     # TODO: convert .iterrows() to .itertuples() for efficiency
#     for i, l in df.iterrows():
#         if l["sample"] in list(analysis_tab["sample"]):
#             outpaths.append(
#                 os.path.join(
#                     outfolder,
#                     "{sample}_{library}.R1.fastq.gz".format(
#                         sample=l["sample"], library=l["library"])
#                 )
#             )
#             outpaths.append(
#                 os.path.join(
#                     outfolder,
#                     "{sample}_{library}.R2.fastq.gz".format(
#                         sample=l["sample"], library=l["library"])
#                 )
#             )
#     return outpaths

def memory_usage_resource() -> float:
    """ Get the current memory usage.

    Returns
    -------
    mem : float
        The memory usage in MB.

    Notes
    -----
    Source: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    """
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # it seems that in OSX the output is different units
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def s_encoding(s: Union[bytes, str]) -> str:
    """ Convert a bytes object to a string, or simply return the given string.

    Parameters
    ----------
    s : Union[bytes, str]
        Input element to convert.

    Returns
    -------
    str
        Input element converted to string (or itself if already string).
    """
    if isinstance(s, bytes):
        return s.decode("utf-8")
    elif isinstance(s, str):
        return s
    else:
        raise TypeError("input not recognised (need to provide"
                        " a str or bytes instance)")


def softclipping(i):
    # TODO: add a proper docstring
    lseq = len(i[9])
    matches = re.findall(r'(\d+)S', i[5])
    sc = sum([int(x) for x in matches])
    return lseq, sc


def get_SAM_header(samfile):
    # is the file compressed?
    if samfile.endswith("gz"):
        samhandle = gzip.GzipFile(samfile, mode='r')
    elif samfile.endswith("bz2"):
        samhandle = bz2.BZ2File(samfile, mode='r')
    else:
        samhandle = open(samfile, 'r')
    comment_count = 0
    header_lines = []
    l = s_encoding(samhandle.readline())
    print(l)
    while l[0] == "@":
        header_lines.append(l)
        comment_count += 1
        l = s_encoding(samhandle.readline())
    return header_lines, comment_count


def check_tmp_dir(fold):
    """ Return the TMP env variable content, or the input folder otherwise. """
    return os.getenv("TMP") or fold


# TODO: re expressions can be simplified
# Functions taken or adapted from assembleMTgenome.py
r = re.compile("#+")
r1 = re.compile("""\^.{1}""")
rr = re.compile("[\+\-]{1}[0-9]+")


# TODO: ref is not used anywhere
def normS(s, ref):
    c = re.finditer(rr, s)
    sl = list(s)
    cc = [(x.start(), x.end()) for x in c]
    for i in cc:
        n = int(''.join(sl[i[0] + 1: i[1]]))
        sl[i[0]: i[1] + n] = ['#' for xx in range(len(sl[i[0]: i[1] + n]))]
    ns = ''.join(sl)
    ns = ns.replace('#', '')
    ss = ''
    for i in ns:
        if i in '.,ACGTNacgtN<>*':
            ss += i
    return ss


def nuc(seq):
    d = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    for i in seq:
        if i in d:
            d[i] += 1
        else:
            d['N'] += 1
    return d


def ff(v, l):
    for i in l:
        x = 0
        for j in i:
            if j in v:
                x += 1
        if x == len(v):
            return i
    return 0


def getIUPAC(f):
    vv = ''.join([i[1] for i in f if i[0] > 0])
    k = ff(vv, DIUPAC.keys())
    if k != 0:
        return DIUPAC[k]
    else:
        return '#'


def freq(d):
    f = []
    for i in d:
        try:
            v = float(d[i]) / sum(d.values())
        except:
            v = 0.0
        f.append((v, i))
    f.sort()
    f.reverse()
    maxv = [f[0]]
    for i in f[1:]:
        if i[0] == maxv[0][0]:
            maxv.append(i)
    if len(maxv) == 1:
        if maxv[0][0] >= CLEV:
            return maxv[0][1]
        else:
            return getIUPAC(f)
    elif len(maxv) > 1:
        return getIUPAC(f)


def get_seq_name(fasta):
    mt_genome = SeqIO.index(fasta, 'fasta')
    if len(mt_genome) != 1:
        sys.exit(("Sorry, but MToolBox at the moment only accepts "
                  "single-contig reference mt genomes."))

    for contig, contig_seq in mt_genome.items():
        seq_name = contig
        seq_length = len(contig_seq)
    return seq_name, seq_length


def nuc_strand(values):
    d = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0}
    for i in values:
        if i[0] in d:
            d[i[0]] += i[1]
        else:
            d['N'] += 1
    return d


# FORMAT CONVERTERS

def pileup2mt_table(pileup=None, ref_fasta=None):
    """
    Returns a list with mt_table_data. To be used to write out the mt-table
    and the fasta output.
    In the MToolBox pipeline, it's used by two rules:
    - pileup2mt_table, which passes the output to write_mt_table function
        to output the mt_table file;
    - make_single_VCF, which passes the output to mt_table_handle2gapped_fasta
        function to generate a gapped fasta of the assembled mt genome.
    """
    # generate mt_table backbone
    # get mt sequence data from ref genome fasta file
    mtdna = {}
    mt_genome = SeqIO.index(ref_fasta, 'fasta')
    if len(mt_genome) != 1:
        sys.exit(("Sorry, but MToolBox at the moment only accepts "
                  "single-contig reference mt genomes."))
    for contig, contig_seq in mt_genome.items():
        for pos, nt in enumerate(contig_seq.seq):
            mtdna[pos + 1] = (nt, ['#', (0, 0, 0, 0), 0, 0.0,
                                   (0, 0, 0, 0, 0, 0, 0, 0)])

    # open input file
    with open(pileup, "r") as f:

        # iterate over pileup
        for i in f:
            if i.strip() == '':
                continue
            l = (i.strip()).split('\t')
            pos = int(l[1])
            if len(l) == 6:
                # TODO: update normS because l[2] is not used anywhere
                ref, seq, qual = l[2], normS(re.sub(r1, "", l[4]), l[2]), l[5]
                s, q = '', 0
                d = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
                     'a': 0, 'c': 0, 'g': 0, 't': 0}
                for j in range(len(seq)):
                    if seq[j] not in '<>*' and ord(qual[j])-33 >= MQUAL:
                        if seq[j] == ".":
                            d[ref.upper()] += 1
                            s += ref.upper()
                        elif seq[j] == ",":
                            d[ref.lower()] += 1
                            s += ref.upper()
                        elif seq[j] in 'acgtACGT':
                            d[seq[j]] += 1
                            s += seq[j].upper()
                        else:
                            pass

                        q += (ord(qual[j]) - 33)
                try:
                    mq = float(q) / len(s)
                except:
                    mq = 0.0
                dnuc = nuc(s)
                mfreq = freq(dnuc)
                lnuc = (dnuc['A'], dnuc['C'], dnuc['G'], dnuc['T'])
                str_nuc = (d['A'], d['C'], d['G'], d['T'],
                           d['a'], d['c'], d['g'], d['t'])
                cnuc = '#'
                if len(s) >= COV:
                    cnuc=mfreq

                mtdna[pos][1][0] = cnuc
                mtdna[pos][1][1] = lnuc
                mtdna[pos][1][2] = len(s)
                mtdna[pos][1][3] = mq
                mtdna[pos][1][4] = str_nuc
            else:
                mtdna[pos][1][0] = '#'

    return mtdna

# End of functions taken from assembleMTgenome
def parse_coverage_data_file(coverage_data_file=None):
    coverage_data = {}
    sam_cov = open(coverage_data_file, 'r')
    for l in sam_cov:
        ref, pos, cov = l.split()
        coverage_data[int(pos)] = int(cov)

    sam_cov.close()
    
    return coverage_data

def sam_cov_handle2gapped_fasta(coverage_data_file=None, ref_mt=None):
    """
    sam_cov_data is a dict {pos : DP}.
    We currently assume ref is only 1 seq.
    """
    ref = SeqIO.index(ref_mt, 'fasta')
    ref_seq = ref[list(ref.keys())[0]].seq
    gapped_fasta = ""
    sam_cov_data = parse_coverage_data_file(coverage_data_file)
    for n in range(len(sam_cov_data)):
        if sam_cov_data[n+1] >= COV:
            gapped_fasta += ref_seq[n]
        else:
            gapped_fasta += "#"
    return gapped_fasta


def write_mt_table(mt_table_data=None, mt_table_file=None):
    """
    Writes out the mt_table file.
    """
    with open(mt_table_file, "w") as mt_table_handle:
        mt_table_handle.write(
            ('Position\tRefNuc\tConsNuc\tCov\tMeanQ\tBaseCount(A,C,G,T)'
             '\tStrandCount(A,C,G,T,a,c,g,t)\n')
        )
        # TODO: these are not used anywhere
        assb, totb = 0, 0
        cop = 0
        maxCval = 1
        for i in range(len(mt_table_data)):
            line = [str(i + 1),
                    mt_table_data[i + 1][0],
                    mt_table_data[i + 1][1][0],
                    str(mt_table_data[i + 1][1][2]),
                    "0:.2f".format(mt_table_data[i + 1][1][3]),
                    str(mt_table_data[i + 1][1][1])]
            line = [str(i + 1),
                    mt_table_data[i + 1][0],
                    mt_table_data[i + 1][1][0],
                    str(mt_table_data[i + 1][1][2]),
                    "0:.2f".format(mt_table_data[i + 1][1][3]),
                    str(mt_table_data[i + 1][1][1]),
                    str(mt_table_data[i + 1][1][4])]
            mt_table_handle.write('\t'.join(line) + '\n')


def mt_table_handle2gapped_fasta(mt_table_data=None):
    """ Generate gapped fasta (string) from mt_table. """
    assb, totb = 0, 0
    cop = 0
    aseq = ''
    for i in range(len(mt_table_data)):
        # if variant is not #, contigs will have reference,
        # otherwise the # that will be subsequently substituted with N
        if mt_table_data[i+1][1][0] != '#':
            aseq += mt_table_data[i+1][0]
        else:
            aseq += mt_table_data[i+1][1][0]
        totb += 1
        # keep this for now
        if mt_table_data[i+1][1][0] != '#':
            assb += 1
            cop += mt_table_data[i+1][1][2]
    return aseq


def gapped_fasta2contigs(gapped_fasta=None):
    """
    Breaks a #-gapped fasta (string) into a contig list, eg:

            gapped_fasta = "ATGCTGTGATTACGTACTG##########CAGTATGTGACGT" -->

        --> contigs = [((1, 19), 'ATGCTGTGATTACGTACTG'), ((30, 42), 'CAGTATGTGACGT')]

    According to a minimum gap length threshold (GLEN). (1, 19) and (30, 42)
    indicate the start and end of contigs on the gapped_fasta
    (then, in turn, on the reference mt genome).
    At the moment, gap length is hardcoded (GLEN = 10).
    """
    # TODO: re expression can be simplified (and imported from above)
    r = re.compile("#+")
    # TODO: r1 and rr are not used anywhere
    r1 = re.compile("""\^.{1}""")
    rr = re.compile("[\+\-]{1}[0-9]+")
    gaps = []
    aseq = gapped_fasta
    fseq = aseq.replace('#', 'N')
    for i in re.finditer(r, gapped_fasta):
        cc = (i.start() + 1, i.end())
        if (cc[1]-cc[0])+1 >= GLEN:
            gaps.append(cc)
    contigs = []
    if len(gaps) != 0:
        for i in range(len(gaps)-1):
            cc = (gaps[i][1]+1, gaps[i+1][0]-1)
            contigs.append((cc, fseq[cc[0]-1:cc[1]]))
        if gaps[0][0] != 1:
            cc = (1, gaps[0][0]-1)
            contigs.insert(0, (cc, fseq[cc[0]-1:cc[1]]))
        if gaps[-1][1] != len(aseq):
            cc = (gaps[-1][1]+1, len(aseq))
            contigs.append((cc, fseq[cc[0]-1:cc[1]]))
        contigs.sort()
    else:
        cc = (1, len(aseq))
        contigs = [(cc, fseq[cc[0]-1:cc[1]+1])]
    return contigs


# TODO: do_softclipping is not used anywhere
def sam_to_fastq(samfile=None, outmt1=None, outmt2=None, outmt=None,
              do_softclipping=True):
    """
    Extract reads from SAM file.
    Reads with softclipping > 1/3 of their length are discarded since could
    be potential non-reference NumtS.
    """
    # count of reads discarded because of softclipping > threshold (hardcoded 1/3)
    sclipped = 0
    print('Extracting FASTQ from SAM...')
    dics = {}
    c = 0

    with gzip.open(outmt, "wb") as mtoutfastq, \
        gzip.open(outmt1, "wb") as mtoutfastq1, \
        gzip.open(outmt2, "wb") as mtoutfastq2, \
        gzip.open(samfile, "rb") as f:
        for i in f:
            i = i.decode("utf-8")
            c += 1
            if c % 100000 == 0:
                print("{} SAM entries processed.".format(c))
            if i.strip() == "" or i.startswith("@"):
                continue
            l = (i.strip()).split("\t")
            if l[2] == "*":
                continue
            lseq, sc = softclipping(l)
            if sc > float(lseq)/3:
                sclipped += 1
                # if soft-clipped read greater than a third of read length,
                # discard the read
                continue
            if len(dics) == 0:
                dics[l[0]] = [l]
            else:
                if l[0] in dics:
                    dics[l[0]].append(l)
                else:
                    # check if reads go to single or paired end file
                    # check if each read in a pair goes to R1 or R2
                    if len(dics) != 1:
                        sys.exit("read data not valid: {}".format(dics))
                    k = [key for key in dics][0]
                    ll = dics[k]
                    if len(ll) == 1:
                        strand, seq, qual = int(ll[0][1]) & 16, ll[0][9], ll[0][10]
                        if strand == 16:
                            seq, qual = reverse_complement(seq), qual[::-1]
                        entry = "\n".join(["@" + ll[0][0], seq, "+", qual]) + "\n"
                        mtoutfastq.write(entry.encode("utf-8"))
                    else:
                        strand, seq, qual = int(ll[0][1]) & 16, ll[0][9], ll[0][10]
                        if strand == 16:
                            seq, qual = reverse_complement(seq), qual[::-1]
                        entry = "\n".join(["@" + ll[0][0], seq, "+", qual]) + "\n"
                        mtoutfastq1.write(entry.encode("utf-8"))
                        strand, seq, qual = int(ll[1][1]) & 16, ll[1][9], ll[1][10]
                        if strand == 16:
                            seq, qual = reverse_complement(seq), qual[::-1]
                        entry = "\n".join(["@" + ll[1][0], seq, "+", qual]) + "\n"
                        mtoutfastq2.write(entry.encode("utf-8"))
                    # create new dics with new read ID
                    dics = {l[0]: [l]}

    return sclipped

#b = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]

def collect_bitwise_flags(n, b=[0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]):
    """
    | Integer  	|    Binary    	|                                Description (Paired Read Interpretation)                               	|
|:--------:	|:------------:	|:-----------------------------------------------------------------------------------------------------:	|
|     1    	|       1      	|                   template having multiple templates in sequencing (read is paired)                   	|
|     2    	|      10      	|          each segment properly aligned according to the aligner (read mapped in proper pair)          	|
|     4    	|      100     	|                                   segment unmapped (read1 unmapped)                                   	|
|     8    	|     1000     	|                         next segment in the template unmapped (read2 unmapped)                        	|
|    16    	|     10000    	|                      SEQ being reverse complemented (read1 reverse complemented)                      	|
|    32    	|    100000    	|    SEQ of the next segment in the template being reverse complemented (read2 reverse complemented)    	|
|    64    	|    1000000   	|                              the first segment in the template (is read1)                             	|
|    128   	|   10000000   	|                              the last segment in the template (is read2)                              	|
|    256   	|   100000000  	|                                         not primary alignment                                         	|
|    512   	|  1000000000  	|                                     alignment fails quality checks                                    	|
|   1024   	|  10000000000 	|                                        PCR or optical duplicate                                       	|
|   2048   	| 100000000000 	| supplementary alignment (e.g. aligner specific, could be a portion of a split read or a tied region)  	|
    """
    flags = []
    for i in b:
        if n & i == 0 and n == 0:
            flags.append(0)
        if n & i != 0:
            flags.append(i)
    return set(flags)
#gzip.open(outmt2, "wb") as mtoutfastq2, \

def sam_to_ids(samfile=None, outmt1=None, outmt=None, keep_orphans=True, return_dict=False, return_files=True):
    """Parses sam file and collect read ids.
    
    Args:
        samfile:            a gzip-compressed SAM file
        outmt:              file with list of SE reads
        outmt1:             file with list of PE reads
        keep_orphans:       wanna keep PE reads who lost their mate?
        return_dict:        wanna return read dict (for debugging)?
        return_files:       wanna return output files (for the pipeline)?
    
    Return:
        if return_dict:     read dict
    """
    # TODO:
    # - set a better log
    # - test!
    c = 0
    f = gzip.open(samfile, "rt") 
    if return_dict:
        read_bitwiseflag_decomp = {}
    if return_files:
        mtoutfastq = gzip.open(outmt, "wt")
        mtoutfastq1 = gzip.open(outmt1, "wt")
    for i in f:
        bitwise_status = True
        c += 1
        if c % 100000 == 0:
            print("{} SAM entries processed.".format(c))
        if i.strip() == "" or i.startswith("@"):
            continue
        l = (i.strip()).split("\t")
        if l[2] == "*":
            continue
        bitwise_flags = collect_bitwise_flags(int(l[1]))
        if 2048 in bitwise_flags: # skip supplementary alignments, we've already met this read
            continue
        elif set([1, 64]).issubset(bitwise_flags): # read paired and first in pair
            paired_status = "PE"
            if return_files:
                mtoutfastq1.write("{}\n".format(l[0]))
        elif 0 in bitwise_flags: # unpaired, mapped
            paired_status = "SE"
            if return_files:
                mtoutfastq.write("{}\n".format(l[0]))
        elif keep_orphans:
            if set([1, 8]).issubset(bitwise_flags): # orphan left from alignment stage
                paired_status = "SE"
                if return_files:
                    mtoutfastq.write("{}\n".format(l[0]))
        else:
            print("Couldn't find assignment for {} with bitwise flag {}".format(l[0], l[1]))
            bitwise_status = False
        if return_dict:
            read_bitwiseflag_decomp[l[0]] = SimpleNamespace(readID=l[0], bitwise_flag=int(l[1]),
                                                        bitwise_decomp=bitwise_flags, bitwise_status=bitwise_status, paired_status=paired_status)
    if return_files:
        mtoutfastq.close()
        mtoutfastq1.close()
    f.close()
    if return_dict:
        return read_bitwiseflag_decomp

def run_seqtk_subset(seqfile=None, id_list=None, outseqfile=None):
    shell("seqtk subseq {seqfile} {id_list} > {outseqfile}")