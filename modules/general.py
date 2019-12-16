#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import bz2
import gzip
import os
import re
import resource
import sys
from typing import Union

from Bio import SeqIO
from Bio.Seq import reverse_complement

from modules.constants import CLEV, COV, DIUPAC, GLEN, MQUAL


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
    lseq = len(i[9])
    sc = re.findall(r'(\d+)S', i[5])
    sc = map(lambda x: int(x), sc)
    sc = sum(sc)
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
    return seq_name


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


def sam_cov_handle2gapped_fasta(sam_cov_data=None, ref_mt=None):
    """
    sam_cov_data is a dict {pos : DP}.
    We currently assume ref is only 1 seq.
    """
    ref = SeqIO.index(ref_mt, 'fasta')
    ref_seq = ref[list(ref.keys())[0]].seq
    gapped_fasta = ""
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
def sam2fastq(samfile=None, outmt1=None, outmt2=None, outmt=None,
              do_softclipping=True):
    """
    Extract reads from SAM file.
    Reads with softclipping > 1/3 of their length are discarded since could
    be potential non-reference NumtS.
    """
    # count of reads discarded because of softclipping > threshold (hardcoded 1/3)
    sclipped = 0
    print('Extracting FASTQ from SAM...')
    mtoutsam = samfile
    mtoutfastq = gzip.GzipFile(outmt, 'wb')
    mtoutfastq1 = gzip.GzipFile(outmt1, 'wb')
    mtoutfastq2 = gzip.GzipFile(outmt2, 'wb')
    f = gzip.GzipFile(mtoutsam, 'rb')
    dics = {}
    c = 0
    for i in f:
        i = i.decode("utf-8")
        c += 1
        if c % 100000 == 0:
            print("{} SAM entries processed.".format(c))
        if i.strip() == '' or i.startswith('@'):
            continue
        l = (i.strip()).split('\t')
        if l[2] == '*':
            continue
        lseq, sc = softclipping(l)
        if sc > float(lseq)/3:
            sclipped += 1
            # if soft-clipped read greater than a third of read length, discard the read
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
                    entry = '\n'.join(['@' + ll[0][0], seq, '+', qual]) + '\n'
                    mtoutfastq.write(entry.encode("utf-8"))
                else:
                    strand, seq, qual = int(ll[0][1]) & 16, ll[0][9], ll[0][10]
                    if strand == 16:
                        seq, qual = reverse_complement(seq), qual[::-1]
                    entry = '\n'.join(['@' + ll[0][0], seq, '+', qual]) + '\n'
                    mtoutfastq1.write(entry.encode("utf-8"))
                    strand, seq, qual = int(ll[1][1]) & 16, ll[1][9], ll[1][10]
                    if strand == 16:
                        seq, qual = reverse_complement(seq), qual[::-1]
                    entry = '\n'.join(['@' + ll[1][0], seq, '+', qual]) + '\n'
                    mtoutfastq2.write(entry.encode("utf-8"))
                # create new dics with new read ID
                dics = {l[0]: [l]}
    mtoutfastq.close()
    mtoutfastq1.close()
    mtoutfastq2.close()
    return sclipped
