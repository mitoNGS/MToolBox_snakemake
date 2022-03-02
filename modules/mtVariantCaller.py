#!/usr/bin/env python

"""
Written by Claudia Calabrese - claudia.calabrese23@gmail.com
    and Domenico Simone - dome.simone@gmail.com
"""

from collections import OrderedDict
import sys
import glob
import gzip
import math
import os
import re
from types import SimpleNamespace
import vcf
from modules.general import parse_coverage_data_file

from Bio.bgzf import BgzfWriter
from Bio import SeqIO
import numpy as np
import scipy as sp
import pandas as pd


def extract_mismatches(seq, qs, len_mism, position_in_read):
    '''extract mismatches from read sequence using MD flag
       
       Parameters
       ----------
       start: int
       end: int
       len_mism: int
       position_in_read: type

       Returns
       -------
       mismatch_seq: str
       mismatch_qs: list
'''
    start = position_in_read - 1
    end = start + len_mism
    mismatch_seq = seq[start:end]
    mismatch_qs = qs[start:end]
    mismatch_qs = list(map(lambda x: ord(x)-33, mismatch_qs))
    return mismatch_seq, mismatch_qs


def check_strand(mate):
    '''check whether read is forward or reverse
        
       Parameters
       ----------
       mate: int
 
       Returns
       -------
       strand: str
    '''
    if mate & 16 == 16:
        strand = '-'
    else:
        strand = '+'
    return strand


def parse_sam_row(row):
    """
    Given a SAM row eg

    0           HWI-ST0866:195:D1J58ACXX:3:1108:15229:52423
    1                                                   163
    2                                           NC_001323.1
    3                                                     1
    4                                                    40
    5                                              47M3I50M
    6                                                     =
    7                                                   136
    8                                                   238
    9     AATTTTATTTTTTAACCTAACTCCCCTACTAAGTGTACCCCCCCTT...
    10    @@@FFFFFHHHHHGHH@GEGHJGIJICEGIGIIIJHHGHGHIGBHI...
    11                                               X2:i:0
    12                                         MD:Z:59T5A31
    13                                          RG:Z:sample
    14                                               NH:i:1
    15                                               HI:i:1
    16                                               NM:i:5
    17                                              SM:i:40
    18                                            XM:Z:100M
    19                                              XO:Z:CU
    20                                              XQ:i:40

    will return:

    md = '59T5A31'
    leftmost = 0
    new_seq = seq as in field 9
    new_qs = qs as in field 10
    strand = 163 (bitwise flag)
    bases = ['59', '5', '31'] (bases as in MD)
    nt = ['T', 'A'] (nt as in MD)
    cigar_bases = ['47', '3', '50']
    cigar_nt = ['M', 'I', 'M']

    All SAM fields of interest are expected to be at specific position because
    they are mandatory (among the first 11), except for MD which is an optional
    one and could be found:
    - in position 12 if the SAM file has not been processed by samtools calmd
    - at the end of the line if the SAM file has been processed by samtools calmd
    """
    assert(isinstance(row, str))
    row = [i.strip() for i in row.split()]
    row[1] = int(row[1])
    row[3] = int(row[3])
    row[4] = int(row[4])
    row[7] = int(row[7])
    row[8] = int(row[8])
    # Find MD field among optional fields
    for field in row[11:]:
        if field.startswith("MD"):
            md = field.split(':')[2]
            if '*' in md:
                sys.stderr.write('SAM field without MD flag or with non-canonical MD flag found. Skip this row\n')
                md = '0'
    leftmost = row[3]-1
    read_id = row[0]
    seq = list(row[9])
    qs = list(row[10])
    strand = check_strand(int(row[1]))
    bases = re.split('[a-zA-Z]', md)
    bases = list(map(lambda x: x.strip('^'), bases))
    bases = filter(None, bases)
    bases = list(map(lambda x: int(x), bases))
    nt = list(filter(None, re.split('[0-9]', md)))
    cigar = row[5]
    cigar_bases = list(filter(None, re.split('[a-zA-Z]', row[5])))
    cigar_nt = list(filter(None, re.split('[0-9]', row[5])))
    new_seq = seq
    new_qs = qs
    return (md, leftmost, new_seq, new_qs, strand, bases, nt, cigar, cigar_bases, cigar_nt)


def read_length_from_cigar(cigar_bases, cigar_nt):
    """ Computes the effective read length of a CIGAR operator taking soft-clipped and deletions into account

    Parameters
    ----------
    cigar_nt: list
              e.g. ['S', 'M', 'D', 'M', 'I', 'M']
    cigar_bases: list
              e.g. [10,  30,   5,  15,  10,  20]
    
    Returns
    -------
    eff_read_length: int
    """
    eff_read_length = 0
    for x, i in enumerate(cigar_nt):
        if i not in ['S', 'D']:
            eff_read_length += int(cigar_bases[x])
    return eff_read_length


def parse_mismatches_from_cigar_md(sam_record, minqs=25, tail=5,
                                   tail_mismatch=5):
    """Extracts mismatch substitutions using MD SAM flag
    
    - MD flag reflects the **mapped portion of the read**  - no soft clipping no
        insertions
    - CIGAR flag reflects the absolute reads length - including soft clipping
        and insertions;
    - The script first equals the length of the read to that of the mapped
        portion and then extracts the variants using the MD flag.
 
    Parameters
    ----------
    sam_record: pandas series
    minqs: int
    tail: int
    tail_mismatch: 5

    Returns
    -------
    list of values
    """
    (md, leftmost, new_seq, new_qs, strand, bases, nt, cigar, cigar_bases, cigar_nt) = parse_sam_row(sam_record)
    # Calculate effective read length (cigar without S and D)
    eff_read_length = read_length_from_cigar(cigar_bases, cigar_nt)
    ins_pos_in_seq = 0
    for n in range(len(cigar_nt)):
        if cigar_nt[n] not in ['I', 'D', 'S', 'H']:
            # if no indel or soft-/hard-clipping, the starting position in the MD is the same
            ins_pos_in_seq += int(cigar_bases[n])
        # cut out from the read sequence and qs the bases representing the ins
        elif cigar_nt[n] == 'I':
        # if insertion found then insert in the read sequence and qs the bases representing the ins
            ins_len = int(cigar_bases[n])
            new_seq = new_seq[:ins_pos_in_seq]+new_seq[(ins_pos_in_seq+ins_len):]
            new_qs = new_qs[:ins_pos_in_seq]+new_qs[(ins_pos_in_seq+ins_len):]
        elif cigar_nt[n] == 'D': #TODO
        # if deletion found then 
            ins_len = int(cigar_bases[n])
            new_seq = new_seq[:ins_pos_in_seq] + ["I"]*ins_len + new_seq[ins_pos_in_seq:]
            new_qs = new_qs[:ins_pos_in_seq] + ["I"]*ins_len + new_qs[ins_pos_in_seq:]
        elif cigar_nt[n] == 'S':
            # if the softclipping is at the beginning of the read
            if n == 0:
                soft_clipped_bases = int(cigar_bases[n])
                new_seq = new_seq[soft_clipped_bases:]
                new_qs = new_qs[soft_clipped_bases:]
            else:
                soft_clipped_bases = int(cigar_bases[n])
                diff = len(new_seq) - soft_clipped_bases
                new_seq = new_seq[0:diff]
                new_qs = new_qs[0:diff]
        else:
            pass
    else:
        pass

    z_pos_evs = list(zip(bases, nt))
    z_pos_evs_ref = []
    z_pos_evs_read = []
    for x, i in enumerate(z_pos_evs):
        if x > 0:
            z_pos_evs_read.append((i[0] + z_pos_evs_read[x-1][0] +
                                   len(z_pos_evs_read[x-1][1].replace('^', '')),
                                   i[1]))
            z_pos_evs_ref.append((i[0] + z_pos_evs_ref[x-1][0] +
                                  len(z_pos_evs_read[x-1][1].replace('^', '')),
                                  i[1]))
        else:
            z_pos_evs_read.append((i[0], i[1]))
            z_pos_evs_ref.append((i[0]+leftmost+1, i[1]))

    z_pos_evs_ref = [i for i in z_pos_evs_ref if "^" not in i[1]]
    z_pos_evs_read = [i for i in z_pos_evs_read if "^" not in i[1]]
    positions_ref = [i[0] for i in z_pos_evs_ref]
    bases_ref = [i[1] for i in z_pos_evs_ref]
    positions_read = [i[0] for i in z_pos_evs_read]
    positions_ref_final = []
    positions_read_final = []
    all_ref = []
    all_mism = []
    all_qs = []
    for x, t in enumerate(positions_read):
        # filter out variants with qs < threshold
        # found some cases where there is a mismatch in soft-clipped zone, this will raise an error
        try:
            if ord(new_qs[t])-33 >= minqs:
                if t >= tail_mismatch and (eff_read_length-t) >= tail_mismatch:
                    positions_ref_final.append(positions_ref[x])
                    positions_read_final.append(positions_read[x])
                    all_ref.append(bases_ref[x])
                    all_mism.append(new_seq[t])
                    all_qs.append(ord(new_qs[t])-33)
        except IndexError: #TODO  - shouldn't we raise here a more human readable error?
            pass
    return positions_ref_final, positions_read_final, all_ref, all_mism, all_qs, strand

def allele_strand_counter(strand):
    """ Initialize a strand counter instance for mismatch detection. 
         
    Parameters
    ----------
    strand: str

    Returns
    -------
    l: list
    """
    if strand == "+":
        l = [1, 0]
    elif strand == "-":
        l = [0, 1]
    return l


def allele_strand_updater(l, allele_strand_count=None):
    """ Updates a strand counter instance for mismatch detection. 
    
    Parameters
    ----------
    l: list
    allele_strand_count: type, default=None
   
    Returns
    -------
    allele_strand_count_new: list
    """
    if allele_strand_count is None:
        allele_strand_count = []
    allele_strand_count_new = []
    for x, i in enumerate(allele_strand_count):
        allele_strand_count_new.append(i + l[x])
    return allele_strand_count_new

def get_per_strand_read_depth(df,rleft, genotype):
    """ Function to calculate per strand read depth used only for indels 
    
    Parameters
    ----------
    df: pandas dataframe
    rleft: int
    genotype: str
   
    Returns
    -------
    sdr: str
    """
    boolean_vector = (df.rleft == rleft) & (df.genotype.astype(str) == genotype)
    strand=df[boolean_vector]['strand'].values[0]
    o = df[boolean_vector]['read_depth_x'].values.tolist()
    if len(o)==1:
        if strand =="+":
            o.append(0) #if there is no rv read supporting
        else:
            o.insert(0,0) #if there is no fwd read supporting
    sdr = str(o[0])+';'+str(o[1])
    return sdr


def varnames(i):
    """ defines global variables for Indels searching
    
    Parameters
    ----------
    i: list
 
    Returns
    -------
    global variables
    """
    CIGAR = i[5]
    readNAME = i[0]
    seq = i[9]
    qs = i[10]
    refposleft = int(i[3]) - 1
    mate = int(i[1])
    strand = check_strand(mate)
    return CIGAR, readNAME, seq, qs, refposleft, strand


# Heteroplasmic fraction quantification
def heteroplasmy(cov, Covbase):
    try:
        if Covbase >= cov:
            Heteroplasmy = float(cov) / float(Covbase)
            het = round(Heteroplasmy, 3)
            return het
        else:
            return 1.0
    except ZeroDivisionError:
        het = 1.0
        return het


# defines mathematical operations
def sum(left):
    s = 0
    for i in left:
        s += int(i)
    return s


def median(l):
    try:
        if len(l) % 2 != 0:
            median = sorted(l)[int(((len(l) + 1) / 2) - 1)]
        else:
            m1 = sorted(l)[int((len(l) / 2) + 1 - 1)]
            m2 = sorted(l)[int((len(l) / 2) - 1)]
            median = (float(m1) + float(m2)) / 2
        return median
    except ZeroDivisionError:
        return 0


def mean(list):
    try:
        s = sum(list)
        m = float(s) / float(len(list))
        return m
    except ZeroDivisionError:
        m = 0
        return m


# defines function for value errors
def error(list):
    try:
        list.remove('')
    except ValueError:
        pass

def qs_context_check(qs, Variant_list, list_of_flanking_bases, tail, Q):
    '''This function checks the median QS of the bases surrounding the Indel. 
    If the Indel is at a distance below tail from read ends or median right or left qs is below the QS threshold, the Indel will be discarded'''
    qsLeft = []
    qsRight = []
    for i in range(len(Variant_list)):
        if list_of_flanking_bases[i] >= tail:
            qsLeft.append(qs[(list_of_flanking_bases[i]-tail):list_of_flanking_bases[i]])
            qsRight.append(qs[list_of_flanking_bases[i]:(list_of_flanking_bases[i]+tail)])
        else:
            qsLeft.append("delete") #number of flanking leftmost bases is below tail
            qsRight.append("delete") #number of flanking leftmost bases is below tail
    qsL=[]
    qsR=[]
    for q in qsLeft:
        if "delete" not in q:
            median_qs_left = median(list(map(lambda x:(ord(x)-33),q))) #calculate the median qs around 5nt leftmost to variant
            if median_qs_left >= Q:
                qsL.append(median_qs_left)
            else:
                qsL.append("delete")
        else:
            qsL.append("delete")
    for q in qsRight:
        if "delete" not in q:
            median_qs_right = median(list(map(lambda x:(ord(x)-33),q))) #calculate the median qs around 5nt rightmost to variant
            if median_qs_right >= Q:
                qsR.append(median_qs_right)
            else:
                qsR.append("delete")
        else:
            qsR.append("delete")
    qs_median = list(map(lambda x:[x[0],x[1]],zip(qsL,qsR))) #list of tuples
    return qs_median

def indels_results(left_tail, right_tail, tail, Indel, var_type, readNAME, strand, indels_flanking,refposleft,qs,Q):
    res = []
    if len(Indel) >= 1 and left_tail >= tail and right_tail >= tail: #check if the first Del and last Del are far more then X nt (tail) from read ends
        res.append([var_type]*len(Indel))
        res.append([readNAME]*len(Indel))
        res.append(strand*len(Indel))
        res.append(refposleft)
        res.append(Indel)
        res.append(qs_context_check(qs,Indel,indels_flanking,tail,Q)) #keep Del
    elif len(Indel) > 1 and left_tail >= tail and right_tail < tail:
        res.append([var_type]*len(Indel))
        res.append([readNAME]*len(Indel))
        res.append(strand*len(Indel))
        res.append(refposlef)
        res.append(Indel)
        res.append(qs_context_check(qs,Indel,indels_flanking,tail,Q))
        res[-1][-1] = ['delete','delete'] #remove the last indel
    elif len(Indel) > 1 and left_tail < tail and right_tail >= tail:
        res.append([var_type]*len(Indel))
        res.append([readNAME]*len(Indel))
        res.append(strand*len(Indel))
        res.append(refposleft)
        res.append(Indel)
        res.append(qs_context_check(qs,Indel,indels_flanking,tail,Q))
        res[-1][0] = ['delete','delete'] #remove the first indel
    else:
        res.append([var_type]*len(Indel))
        res.append([readNAME]*len(Indel))
        res.append(strand*len(Indel))
        res.append(refposleft)
        res.append(Indel)
        res.append([("delete","delete")]) #remove indel
    res_final =  list(map(lambda x:[x[0],x[1],x[2],x[3],x[4],x[5]],zip(res[0],res[1],res[2],res[3],res[4],res[5]))) #create list of lists of Indels
    #returns a list of values:
    #[SRR043366.13710149', '-', 309, ['T', 'C'], [33, 33]] for Insertions
    #['SRR043366.15373156', '-', 16188, range(16189, 16190), [34, 34]] for Deletions
    return res_final

def get_Final_dictionary(Final, df, vartype):
    '''Final is a dictionary
       df is a pandas dataframe
       vartype == "ins" | "del" '''
    for i,x in df.iterrows():
        srd = get_per_strand_read_depth(df,x.rleft,x.genotype)
        if x.rleft not in Final:
            Final[x.rleft] = [[vartype,x.genotype,x.mean_qs,x.read_depth_y,srd]] #read_depth_y is the total read depth
        else:
            list_genotypes = list(map(lambda x:x[1], Final[x.rleft]))
            n = np.in1d(list_genotypes, x.genotype)
            if sum(n) == 1: #if the genotype is already there
                Final[x.rleft][np.where(n == True)[0][0]][3]+x.read_depth_y #calculate final read depth for that position
            else: #the genotype is not there
                Final[x.rleft].append([vartype,x.genotype,x.mean_qs,x.read_depth_y,srd])
    return Final

# defines the function searching for and filtering indels within the read sequence
def SearchINDELsintoSAM(readNAME,strand,CIGAR,seq,qs,refposleft,tail=5,Q=25): #TODO - change tail and Q to customizable values
    m=re.compile(r'[a-z]', re.I)
    res = []
    #take indexes of operators in CIGAR
    op_start = [x.start() for x in m.finditer(CIGAR)]
    CIGAR_sp = np.array(list(CIGAR))
    all_changes = CIGAR_sp[op_start]
    op_start = np.array(op_start)
    list_of_indexes = [[0,op_start[0]]]
    i = 0
    while i < len(op_start)-1:
        if i == len(op_start)-2:
            t = [op_start[i]+1,op_start[-1]]
            list_of_indexes.append(t)
        else:
            t = [op_start[i]+1,op_start[i+1]]
            list_of_indexes.append(t)
        i += 1
    #slice CIGAR based on start:end in list_of_indexes
    all_bp = np.array(list(map(lambda x:int(CIGAR[x[0]:x[1]]),list_of_indexes)))
    if 'D' in CIGAR or 'N' in CIGAR: #GMAP can use also N for large deletions
        #DELETIONS
        #boolean vector indicating position of D
        bv_del = np.in1d(all_changes,'D') | np.in1d(all_changes,'N')
        var_type = 'Del'
        #boolean vector indicating position of Hard clipped (H) and Soft clipped bases (S) to be removed from leftmost count
        bv_hard_or_soft = (np.in1d(all_changes,'H')) | (np.in1d(all_changes,'S'))
        #dummy vector
        d = np.zeros(len(bv_del))
        #adding leftmost positions, excluding those preceding H and S
        d[~bv_hard_or_soft] = all_bp[~bv_hard_or_soft]
        #calculate cumulative number of bp before each del
        cum_left = np.cumsum(d)
        dels_indexes = np.where((all_changes=='D') | (all_changes=='N'))[0]
        flanking_dels_indexes = dels_indexes-1
        #calculate leftmost positions to dels within the read
        refposleft_dels = cum_left[flanking_dels_indexes]
        refposleft_dels = refposleft_dels + refposleft
        refposleft_dels = refposleft_dels.astype(int).tolist()
        #get Deletion coordinates
        dels_indexes = op_start[bv_del]-1
        dels = list(map(lambda x:int(x),CIGAR_sp[dels_indexes]))
        list_dels = zip(refposleft_dels,dels)
        Del = list(map(lambda x:range(x[0]+1,x[0]+1+x[1]),list_dels))
        #get left and right tails of dels
        dels_flanking = all_bp[flanking_dels_indexes]
        left_tail = dels_flanking[0]
        right_tail = len(seq)-sum(dels_flanking)
        res_del = indels_results(left_tail, right_tail, tail, Del, var_type, readNAME, strand, dels_flanking,refposleft_dels,qs,Q)
        res.extend(res_del)
    if 'I' in CIGAR:
        #INSERTIONS
        ins_indexes = np.where(all_changes=='I')[0]
        bv_ins = np.in1d(all_changes,'I')
        var_type = 'Ins'
        #boolean vector indicating position of Hard clipped (H) and Soft clipped bases (S) to be removed from leftmost count
        bv_hard_or_soft = (np.in1d(all_changes,'H')) | (np.in1d(all_changes,'S'))
        #dummy vector with same length as many ins in the CIGAR
        i = np.zeros(len(bv_ins))
        #adding leftmost positions, excluding those preceding H and S
        i[~bv_hard_or_soft] = all_bp[~bv_hard_or_soft]
        #calculate cumulative number of bp before each ins and getting the flanking index in the read
        i[bv_ins] = 0
        cum_left = np.cumsum(i)
        flanking_ins_indexes = ins_indexes-1
        #calculate leftmost positions to ins within the read
        refposleft_ins = ((cum_left[flanking_ins_indexes])+refposleft).astype(int)
        #get Insertion length
        op_bases = op_start[bv_ins]-1
        ins_length = list(map(lambda x:int(x),CIGAR_sp[op_bases]))
        list_ins = zip(refposleft_ins,ins_length)
        Ins = list(map(lambda x:range(x[0]+1,x[0]+1+x[1]),list_ins))
        ins_flanking = all_bp[flanking_ins_indexes]
        left_tail = ins_flanking[0]
        right_tail = len(seq)-sum(ins_flanking)
        res_ins = indels_results(left_tail, right_tail, 5, Ins, var_type, readNAME, strand, ins_flanking,refposleft_ins,qs,30)
        #this returns a list of lists such as: [['Ins', 'SRR043366.13710149', '-', 309, range(310, 312), [33, 33]]]
        #get quality per Ins using Ins positions relative to the read
        qs = np.array(list(qs))
        seq = np.array(list(seq))
        i = np.zeros(len(bv_ins))
        i[~bv_hard_or_soft] = all_bp[~bv_hard_or_soft]
        if "bv_del" in locals():
            i[bv_del] = 0
        cum_ins = (np.cumsum(i)[ins_indexes]).astype(int) #number of bases before the "I" in CIGAR string
        #ins_cum_bases = cum_ins[ins_indexes].astype(int) DELETE THIS IS WRONG
        #ins_start_position = np.cumsum(i)[0] DELETE THIS
        ins_start_position = cum_ins - ins_length #ins start position relative to read (considering 0-base counting)
        list_ins = zip(ins_start_position,ins_length)
        Ins2 = list(map(lambda x:range(x[0],x[0]+x[1]),list_ins))
        qsInsASCI = list(map(lambda x: qs[x].tolist(),Ins2))
        Ins = list(map(lambda x: ''.join(seq[x].tolist()),Ins2))
        #add to results a list with quality scores of ins
        for x in range(len(qsInsASCI)):
            res_ins[x].append(list(map(lambda x:ord(x)-33,qsInsASCI[x]))) #adding an extra value to the insertion res list with QS of each insertion
            res_ins[x][4] = Ins[x]
        res.extend(res_ins)
    return res

def parse_indels(df, Q, minrd, tag):
    boolean_vector1 = df[tag].astype(str).str.contains('delete')
    df = df[~boolean_vector1]
    #filters on qs
    boolean_vector2 = np.array(list(map(lambda x:x[0]<Q,df[tag])) or list(map(lambda x:x[1]<Q,df[tag])))
    df = df[~boolean_vector2]
    df.insert(5,'mean_qs',list(map(lambda x:np.mean(x),df[tag]))) 
    df['genotype'] = df['genotype'].astype(str)
    df_counts = df.groupby(['rleft','strand', 'genotype']).count().reset_index()[['rleft','strand','genotype','Type']]
    df_counts.columns = ['rleft','strand','genotype','read_depth']
    df_median_qs = df.groupby(['rleft','genotype']).median().reset_index()
    df_median_qs.columns = ['rleft','genotype','mean_qs']
    df_final = pd.merge(df_counts,df_median_qs, on=['rleft','genotype'],how='inner')
    depth = df_counts.groupby(['rleft','genotype'])['read_depth'].sum().reset_index()
    depth = depth[['rleft','genotype','read_depth']]
    rleft_pos_to_keep = np.array(depth[depth['read_depth']>=minrd]['rleft'])
    boolean_vector3 = np.in1d(df_final['rleft'],rleft_pos_to_keep)
    genotype_to_keep = np.array(depth[depth['read_depth']>=minrd]['genotype'])
    boolean_vector4 = np.in1d(df_final['genotype'],genotype_to_keep)
    df_final = df_final[(boolean_vector3 & boolean_vector4)]
    df_final = pd.merge(df_final,depth, on=['rleft','genotype'],how='inner')
    return df_final



# Wilson confidence interval lower bound
def CIW_LOW(het, totrd):
    """ The function calculates the heteroplasmic fraction and the related
        confidence interval with 95% of coverage probability,
        considering a Wilson score interval when n<=40
        CIw = [1/(1+(1/n)*z^2)] * [p + (1/2n)*z^2 +- z(1/n *(p*q) + ((1/(4n^2))*z^2))^1/2]
    """
    p = het
    n = totrd
    z = 1.96
    q = 1-het
    num = p * q
    squarez = z * z
    squaren = n * n
    wilsonci_low = round((p + (z * z) / (2 * n) - z *
                          (math.sqrt(p * q / n + (z * z) / (4 * (n * n))))) /
                         (1 + z * z / n), 3)
    if wilsonci_low < 0.0:
        return 0.0
    else:
        return wilsonci_low


# Wilson confidence interval upper bound
def CIW_UP(het, totrd):
    """ The function calculates the heteroplasmic fraction and the related
        confidence interval with 95% of coverage probability,
        considering a Wilson score interval when n<=40
    CIw = [1/(1+(1/n)*z^2)] * [p + (1/2n)*z^2 +- z(1/n *(p*q) + ((1/(4n^2))*z^2))^1/2]
    """
    p = het
    n = totrd
    z = 1.96
    q = 1-het
    num = p*q
    squarez = z * z
    squaren = n * n
    wilsonci_up = round((p + (z * z) / (2 * n) + z *
                         (math.sqrt(p * q / n + (z * z) / (4 * (n * n))))) /
                        (1 + z * z / n), 3)
    if wilsonci_up > 1.0:
        return 1.0
    else:
        return wilsonci_up


# Agresti-Coull confidence interval lower bound
def CIAC_LOW(rd, totrd):
    """ The function calculates the heteroplasmic fraction and the related
        confidence interval for heteroplasmic fraction with 95% of coverage
        probability, considering the Agresti-Coull interval when n>40.
    """
    z = 1.96
    if rd > totrd:
        totrd = rd
    X = rd + (z*z) / 2
    N = totrd + (z*z)
    P = X / N
    Q = 1 - P
    agresticoull_low = round(P - (z * (math.sqrt(P * Q / N))), 3)
    if agresticoull_low < 0.0:
        return 0.0
    else:
        return agresticoull_low


# Agresti-Coull confidence interval upper bound
def CIAC_UP(rd, totrd):
    """ The function calculates the heteroplasmic fraction and the related
        confidence interval for heteroplasmic fraction with 95% of coverage
        probability, considering the Agresti-Coull interval when n>40.
    """
    z = 1.96
    if rd > totrd:
        totrd = rd
    X = rd + (z*z) / float(2)
    N = totrd + (z*z)
    P = X / N
    Q = 1 - P
    agresticoull_up = round(P + (z * (math.sqrt(P * Q / N))), 3)
    if agresticoull_up > 1.0:
        return 1.0
    else:
        return agresticoull_up


# IUPAC dictionary
dIUPAC = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
          'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'],
          'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
          'N': ['A', 'C', 'G', 'T']}


# searches for IUPAC codes and returns the ambiguity
# returns '' if nucleotide in reference is N
def getIUPAC(ref_var, dIUPAC):
    iupac_code = ['']
    for i in dIUPAC.items():
        i[1].sort()
        if ref_var == i[1]:
            iupac_code = [i[0]]
    return iupac_code


# TODO: duplicated function
# function copied from Snakefile
def s_encoding(s):
    if type(s) == bytes:
        return s.decode("utf-8")
    elif type(s) == str:
        return s

def mtvcf_main_analysis(mtable_file=None, coverage_data_file=None, sam_file=None,
                        name2=None, tail=5, Q=25, minrd=5, ref_mt=None,
                        tail_mismatch=5):

    coverage_data = parse_coverage_data_file(coverage_data_file)

    if sam_file.endswith("gz"):
        sam = gzip.GzipFile(sam_file, mode = 'r')
    else:
        sam = open(sam_file, 'r')
    # TODO: these are not used anywhere
    CIGAR = ''
    readNAME = ''
    seq = ''
    qs = ''
    refposleft = ''
    mate = ''
    # populate:
    # - mtDNA: a list of bases in the reference mt genome
    # - coverage: a list of DP as calculated by samtools depth
    mtDNA = []
    Coverage = []
    ref = SeqIO.index(ref_mt, 'fasta')
    ref_seq = ref[list(ref.keys())[0]].seq
    for n in range(len(ref_seq)):
        Coverage.append(coverage_data[n + 1])
        mtDNA.append(ref_seq[n])
    mtDNAseq = "".join(mtDNA)
    ## apply functions to sam file and write outputs into a dictionary
    # add indels 
    dic = {}
    dic['Ins'] = []
    dic['Del'] = []
    print("\nsearching for indels in {0}.. please wait...\n".format(name2))
    for i in sam:
        i = s_encoding(i)
        if i.startswith("@"):
            continue
        i = i.split('\t')
        [CIGAR, readNAME, seq, qs, refposleft, strand] = varnames(i)
        mm = 0
        if 'I' in CIGAR or 'D' in CIGAR:
            r = SearchINDELsintoSAM(readNAME, strand, CIGAR, seq, qs, refposleft,
                                    tail=tail, Q=Q)
            # r is: [['Ins' or 'Del', readNAME, strand, rLeft, Del/Ins, qsDel/qsIns]] where
            # rleft is the leftmost position to Indel
            # Del/Ins is a list of two values like [int,int] 
            ##WARNING - Del/Ins values have different meaning:
            # if indel is deletion then Del is a range(int,int) which returns a tuple with absolute positins of the deleted nt in the mapped read
            # if indel is insertion then Ins is a list of two values like [str,str]  length consistent with Del
            #############
            # qsDel/qsIns is a list like [int,int] with values indicating the median quality scores of surrounding bases to the indel
            for indel in r:
                dic[indel[0]].append(indel[1:])
    #the following code is used to exclude indels based on sequencing context
    Ins_dict = {'Ins': dic['Ins']}
    Del_dict = {'Del': dic['Del']}
    Ins_list = [[key] + i for key,value in Ins_dict.items() for i in value]
    Del_list = [[key] + i for key,value in Del_dict.items() for i in value]
    df_Ins = pd.DataFrame(Ins_list, columns = ['Type', 'readName', 'strand', 'rleft', 'genotype', 'quality', 'qs_flanking'])
    df_Del = pd.DataFrame(Del_list, columns = ['Type', 'readName', 'strand', 'rleft', 'genotype', 'qs_flanking'])
    df_Ins_final = parse_indels(df_Ins, Q, minrd, 'quality')
    df_Del_final = parse_indels(df_Del, Q, minrd, 'qs_flanking')
    Final = {}           
    Final = get_Final_dictionary(Final,df_Ins_final,'ins')
    Final = get_Final_dictionary(Final,df_Del_final,'del')    
    Indels = {}
    Indels[name2] = []
    for i in Final:
        if len(Final.get(i)) > 0:
            for x in Final.get(i):
                if x[0] == 'ins' and x[1] != []:  # is not empty
                    bases = x[1]
                    qs = round(float(x[2]),2)
                    rd = int(x[3])
                    strand = x[4]
                    Refbase = mtDNAseq[int(i)-1]
                    Variant = Refbase + bases
                    totrd = int(Coverage[int(i)-1])
                    if rd > totrd:
                        sys.stderr.write("insertion in pos {} with per base rd > total rd. Assuming total rd is equal to the bigger value\n".format(str(i)))
                        totrd = rd
                    hetfreq = (heteroplasmy(rd,totrd))
                    if totrd <= 40:
                        het_ci_low = CIW_LOW(hetfreq, totrd)
                        het_ci_up = CIW_UP(hetfreq, totrd)
                    else:
                        het_ci_low = CIAC_LOW(rd, totrd)
                        het_ci_up = CIAC_UP(rd, totrd) 
                    ins = [i, Refbase, totrd, [Variant], [rd], [[strand]], [qs],
                           [hetfreq], [het_ci_low], [het_ci_up], 'ins']
                    Indels[name2].append(ins)
                else:
                    if x[1] != [] :  # is not empty
                        Refbase = []
                        qs = x[2]
                        rd = int(x[3])
                        strand = x[4]
                        deletions = []
                        dels = eval(x[1])
                        delflank = dels[0]-2
                        delfinal = dels[-1]
                        covlist = Coverage[delflank:delfinal]
                        convert = list(map(lambda x: int(x), covlist))
                        totrd = round(np.median(convert),0) #median read depth of the region encompassing the del (samtools)
                        if rd > totrd:
                            sys.stderr.write("deletion in pos {} with per base rd > total rd. Assuming total rd is equal to the bigger value\n".format(str(i)))
                        hetfreq = heteroplasmy(rd, totrd)
                        if totrd <= 40:
                            het_ci_low = CIW_LOW(hetfreq, totrd)
                            het_ci_up = CIW_UP(hetfreq, totrd)
                        else:
                            het_ci_low = CIAC_LOW(rd, totrd)
                            het_ci_up = CIAC_UP(rd, totrd)
                        deletions.append(mtDNAseq[delflank])
                        Refbase.append(mtDNAseq[delflank:delfinal])
                        dele = [(dels[0]-1), Refbase, totrd, deletions, [rd],
                                [[strand]], [qs], [hetfreq], [het_ci_low], [het_ci_up], 'del']
                        Indels[name2].append(dele)
    # Mismatch detection
    print("\n\nsearching for mismatches in {0}.. please wait...\n\n".format(name2))

    mismatch_dict = mismatch_detection(sam=sam_file, coverage_data=coverage_data,
                                       tail_mismatch=tail_mismatch)
    x = 0  # alignment counter
    print("mismatch_dict length before filtering for allele_DP: {}".format(len(mismatch_dict)))
    print(mismatch_dict[j] for j in list(mismatch_dict.keys())[:5])
    for POS in mismatch_dict:
        good_alleles_index = [mismatch_dict[POS].allele_DP.index(i)
                              for i in mismatch_dict[POS].allele_DP if i > 5]
        mismatch_dict[POS].alleles = [mismatch_dict[POS].alleles[j]
                                      for j in good_alleles_index]
        # if this filters out all alleles, delete key from dict
        # if len(mismatch_dict[POS].alleles) > 0:
        mismatch_dict[POS].allele_DP = [mismatch_dict[POS].allele_DP[j]
                                        for j in good_alleles_index]
        mismatch_dict[POS].allele_strand_count = [mismatch_dict[POS].allele_strand_count[j]
                                                  for j in good_alleles_index]
    mismatch_dict = {POS: data for POS, data in mismatch_dict.items()
                     if len(data.alleles) > 0}
    print("mismatch_dict length after filtering for allele_DP: {}".format(len(mismatch_dict)))
    print(mismatch_dict[j] for j in list(mismatch_dict.keys())[:5])
    # for now the mismatch dict goes into the Subst dict
    Subst = {}
    Subst[name2] = []
    for POS in mismatch_dict:
        mismatch_dict[POS].hetfreq = list(map(lambda x: heteroplasmy(x,
                                                                     mismatch_dict[POS].DP),
                                              mismatch_dict[POS].allele_DP))
        if mismatch_dict[POS].DP <= 40:
            mismatch_dict[POS].het_ci_low = list(map(lambda x: CIW_LOW(x,
                                                                       mismatch_dict[POS].DP),
                                                     mismatch_dict[POS].hetfreq))
            mismatch_dict[POS].het_ci_up = list(map(lambda x: CIW_UP(x, mismatch_dict[POS].DP),
                                                    mismatch_dict[POS].hetfreq))
        else:
            mismatch_dict[POS].het_ci_low = list(map(lambda x: CIAC_LOW(x,
                                                                        mismatch_dict[POS].DP),
                                                     mismatch_dict[POS].allele_DP))
            mismatch_dict[POS].het_ci_up = list(map(lambda x: CIAC_UP(x,
                                                                      mismatch_dict[POS].DP),
                                                    mismatch_dict[POS].allele_DP))
        a = [POS, mismatch_dict[POS].REF, mismatch_dict[POS].DP,
             mismatch_dict[POS].alleles, mismatch_dict[POS].allele_DP,
             join_allele_strand_count(mismatch_dict[POS].allele_strand_count), 'PASS',
             mismatch_dict[POS].hetfreq, mismatch_dict[POS].het_ci_low,
             mismatch_dict[POS].het_ci_up, 'mism']
        Subst[name2].append(a)

    Indels[name2].extend(Subst[name2])
    return Indels  # it's a dictionary

def join_allele_strand_count(allele_strand_count):
    """A patch to get allele strand counts for mismatches in the same format as those for indels."""
    f = []
    for j in allele_strand_count:
        f.append([";".join([str(i) for i in j])])
    return f

def mismatch_detection(sam=None, coverage_data=None, tail_mismatch=5):
    if sam.endswith("gz"):
        sam_handle = gzip.GzipFile(sam, mode='r')
    else:
        sam_handle = open(sam, 'r')

    mismatch_dict = {}
    for r in sam_handle:
        if r.startswith("@"):
            continue
        (positions_ref, positions_read, all_ref, all_mism,
         all_qs, strand) = parse_mismatches_from_cigar_md(r,
                                                          tail_mismatch=tail_mismatch)
        if positions_ref == []:
            continue
        for mut in zip(positions_ref, positions_read, all_ref, all_mism, all_qs):
            POS = mut[0]
            REF = mut[2]
            allele = mut[3]
            if POS in mismatch_dict:
                try:
                    # check if that allele has already been found for that position
                    allele_index = mismatch_dict[POS].alleles.index(allele)
                    mismatch_dict[POS].allele_DP[allele_index] += 1
                    mismatch_dict[POS].allele_strand_count[allele_index] = allele_strand_updater(
                        l=allele_strand_counter(strand),
                        allele_strand_count=mismatch_dict[POS].allele_strand_count[allele_index]
                    )
                except:
                    allele_index = len(mismatch_dict[POS].alleles)
                    mismatch_dict[POS].alleles.append(allele)
                    mismatch_dict[POS].allele_DP.append(1)
                    mismatch_dict[POS].allele_strand_count.append(allele_strand_counter(strand))
            else:
                # DP needs to be parsed from bcftools/bedtools output
                mismatch_dict[POS] = SimpleNamespace(POS=POS, REF=REF,
                                                     DP=coverage_data[POS],
                                                     alleles=[allele],
                                                     allele_DP=[1],
                                                     allele_strand_count=[allele_strand_counter(strand)])
    sam_handle.close()
    return mismatch_dict

def get_consensus_single(i, hf_max=0.8, hf_min=0.2):
    consensus_value = []
    if len(i) != 0:
        for var in i:
            if var[-1] == 'mism' and max(var[7]) > hf_max:
                index = var[7].index(max(var[7]))
                basevar = var[3][index]
                res = [var[0], [basevar], 'mism']
                consensus_value.append(res)
            elif var[-1] == 'mism' and max(var[7]) >= hf_min and max(var[7]) <= hf_max:
                basevar=np.array([var[1]] + var[3])
                # keep only basevar >= hf_min for IUPAC
                ref_hf = 1-np.sum(var[7])
                hf_var = [ref_hf]
                hf_var.extend(var[7])
                hf_var = np.array(hf_var)
                ii = np.where(hf_var >= hf_min)[0]
                basevar = basevar[ii].tolist()
                basevar.sort()
                a = getIUPAC(basevar, dIUPAC)
                res = [var[0], a, 'mism']
                consensus_value.append(res)
            elif var[-1] == 'mism' and max(var[7]) < hf_min:  # put the reference allele in consensus
                res = [var[0], [var[1]], 'mism']
            elif var[-1] == 'ins' and max(var[7]) > hf_max:
                index = var[7].index(max(var[7]))
                basevar = var[3][index]
                res = [var[0], [basevar], 'ins']
                consensus_value.append(res)
            elif var[-1] == 'del' and max(var[7]) > hf_max:
                index = var[7].index(max(var[7]))
                basevar = var[3][index]
                del_length = len(var[1][0]) - len(basevar)
                start_del = var[0] + 1
                end_del = start_del + del_length
                res = [var[0], range(start_del, end_del), 'del']
                consensus_value.append(res)
            else:
                pass
        return consensus_value


def get_consensus(dict_of_dicts, hf_max, hf_min):
    """ Dictionary of consensus variants, for fasta sequences. """
    Consensus = {}
    for i in dict_of_dicts:
        Consensus[i] = get_consensus_single(dict_of_dicts[i], hf_max, hf_min)
    return Consensus


# TODO: lots of type comparison to replace with isinstance
def VCFoutput(dict_of_dicts, reference='mt_genome', vcffile='sample',
              seq_name='seq', seq_length=0):
    print("Reference sequence used for VCF: {}".format(reference))
    print("Seq_name is {}".format(seq_name))
    VCF_RECORDS = []
    present_pos = set()
    # for each sample in dict_of_dicts
    for sample in dict_of_dicts.keys():
        # gets variants found per sample
        val = dict_of_dicts[sample]
        c = 0
        for variant in val:
            # variant[5] is SDP
            # if the v. position was never encountered before, is heteroplasmic and is a deletion
            if variant[0] not in present_pos and max(variant[7]) < 1 and variant[-1] == 'del':
                allelecount = [1] * len(variant[1])
                aplotypes = list(map(lambda x: x+1, range(len(allelecount))))
                r = vcf.model._Record(CHROM=seq_name, POS=variant[0], ID='.',
                                      REF=variant[1], ALT=variant[3], QUAL='.',
                                      FILTER='PASS',
                                      INFO=OrderedDict([('AC', allelecount),
                                                        ('AN', len(variant[1])+1)]),
                                      FORMAT='GT:DP:HF:CILOW:CIUP:SDP',
                                      sample_indexes={sample: ''}, samples=[])
                r._sample_indexes[sample] = [[0] + aplotypes, variant[2],
                                             variant[7], variant[8], variant[9],
                                             variant[5]]
                r.samples.append(sample)
                if len(variant[3]) > 1:
                    r.REF = r.REF * len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR = [variant[-1]] * len(variant[3])
            # if the v. position was never encountered before, is heteroplasmic and is not a deletion
            elif variant[0] not in present_pos and max(variant[7]) < 1 and variant[-1] != 'del':
                allelecount = [1] * len(variant[3])
                aplotypes = list(map(lambda x: x+1, range(len(allelecount))))
                r = vcf.model._Record(CHROM=seq_name, POS=variant[0], ID='.',
                                      REF=[variant[1]], ALT=variant[3], QUAL='.',
                                      FILTER='PASS',
                                      INFO=OrderedDict([('AC', allelecount),
                                                        ('AN', len(variant[3])+1)]),
                                      FORMAT='GT:DP:HF:CILOW:CIUP:SDP',
                                      sample_indexes={sample: ''}, samples=[])
                r._sample_indexes[sample] = [[0] + aplotypes, variant[2],
                                             variant[7], variant[8], variant[9],
                                             variant[5]]
                r.samples.append(sample)
                if len(variant[3]) > 1:
                    r.REF = r.REF * len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR = [variant[-1]] * len(variant[3])
            # if the v. position was never encountered before,is homoplasmic and is a deletion
            elif variant[0] not in present_pos and max(variant[7]) >= 1 and variant[-1] == 'del':
                allelecount = [1] * len(variant[1])
                r = vcf.model._Record(CHROM=seq_name, POS=variant[0], ID='.',
                                      REF=variant[1], ALT=variant[3], QUAL='.',
                                      FILTER='PASS',
                                      INFO=OrderedDict([('AC', allelecount),
                                                        ('AN', 1)]),
                                      FORMAT='GT:DP:HF:CILOW:CIUP:SDP',
                                      sample_indexes={sample: ''}, samples=[])
                r._sample_indexes[sample] = [[1], variant[2], variant[7],
                                             variant[8], variant[9], variant[5]]
                r.samples.append(sample)
                if len(variant[3])>1:
                    r.REF = r.REF * len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR = [variant[-1]] * len(variant[3])
            # if the v. position was never encountered before,is homoplasmic and is not a deletion
            elif variant[0] not in present_pos and max(variant[7]) >= 1 and variant[-1] != 'del':
                allelecount = [1] * len(variant[3])
                r = vcf.model._Record(CHROM=seq_name, POS=variant[0], ID='.',
                                      REF=[variant[1]], ALT=variant[3], QUAL='.',
                                      FILTER='PASS',
                                      INFO=OrderedDict([('AC', allelecount),
                                                        ('AN', 1)]),
                                      FORMAT='GT:DP:HF:CILOW:CIUP:SDP',
                                      sample_indexes={sample: ''}, samples=[])
                r._sample_indexes[sample] = [[1], variant[2], variant[7],
                                             variant[8],variant[9], variant[5]]
                r.samples.append(sample)
                if len(variant[3]) > 1:
                    r.REF = r.REF * len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR = [variant[-1]] * len(variant[3])
            # If the v.position was encountered before
            elif variant[0] in present_pos and max(variant[7]) < 1:
                for i in VCF_RECORDS:
                    if variant[0] == i.POS:
                        # when there are multiple variants for a position of the same individual
                        if sample in i.samples and type(variant[1]) == type(list()):
                            for x in range(len(variant[3])):
                                if variant[3][x] in i.ALT and variant[1][x] in i.REF:
                                    index = i.ALT.index(variant[3][x])
                                    i.INFO['AC'][index] += 1
                                    i.INFO['AN'] += 1
                                    aplotype = index + 1
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                elif variant[3][x] in i.ALT and variant[1][x] not in i.REF:
                                    i.INFO['AC'].append(1)
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                    i.INFO['AN'] += 1
                                    index = len(i.ALT)-1  # the alt allele added to i.ALT is the last index
                                    aplotype = len(i.INFO['AC'])
                                    i.TYPEVAR.append(variant[-1])
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                else:
                                    i.REF.append(variant[1][x])
                                    i.ALT.append(variant[3][x])
                                    i.INFO['AC'].append(1)
                                    i.INFO['AN'] += 1
                                    index = i.ALT.index(variant[3][x])
                                    i.TYPEVAR.append(variant[-1])
                                    aplotype = index + 1
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                        # for multiple variants of a position in different individuals
                        elif sample not in i.samples and type(variant[1]) == type(list()):
                            i.INFO['AN'] += 1
                            i.samples.append(sample)
                            for x in range(len(variant[3])):
                                if variant[3][x] in i.ALT and variant[1][x] in i.REF:
                                    index = i.REF.index(variant[1][x])
                                    i.INFO['AC'][index] += 1
                                    i.INFO['AN'] += 1
                                    aplotype = index + 1
                                    genotype = [aplotype]
                                elif variant[3][x] in i.ALT and variant[1][x] not in i.REF:
                                    i.INFO['AC'].append(1)
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                    i.INFO['AN'] += 1
                                    index = i.REF.index(variant[1][x])
                                    aplotype = len(i.INFO['AC'])
                                    genotype = [aplotype]
                                    i.TYPEVAR.append(variant[-1])
                                else:
                                    i.INFO['AC'].append(1)
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                    i.INFO['AN'] += 1
                                    index = i.ALT.index(variant[3][x])
                                    aplotype = index + 1
                                    genotype = [aplotype]
                                    i.TYPEVAR.append(variant[-1])
                            i._sample_indexes.setdefault(sample,
                                                         [[0] + genotype,
                                                          variant[2], variant[7],
                                                          variant[8], variant[9],
                                                          variant[5]])
                        elif sample in i.samples and type(variant[1]) != type(list()):
                            for allele in variant[3]:
                                if allele not in i.ALT:
                                    i.REF.append(variant[1])
                                    i.ALT.append(allele)
                                    i.INFO['AC'].append(1)
                                    i.INFO['AN'] += 1
                                    index = i.ALT.index(allele)
                                    aplotype = index + 1
                                    hf_index = variant[3].index(allele)
                                    if type(i._sample_indexes[sample][0]) == type(list()):
                                        i._sample_indexes[sample][0].append(aplotype)
                                    else:
                                        hap = i._sample_indexes[sample][0]
                                        i._sample_indexes[sample][0] = [hap].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][hf_index])
                                    i._sample_indexes[sample][3].append(variant[8][hf_index])
                                    i._sample_indexes[sample][4].append(variant[9][hf_index])
                                    if i.POS == 3107 and variant[10] == 'mism':
                                        # this is specific for human genome
                                        # not hardcodable, maybe have some special way to get this info
                                        i._sample_indexes[sample][5].append(variant[5][hf_index])
                                    else:
                                        i._sample_indexes[sample][5].append(variant[5][hf_index])
                                    i.TYPEVAR.append(variant[-1])
                                else:
                                    index = i.ALT.index(allele)
                                    i.INFO['AC'][index] += 1
                                    i.INFO['AN'] += 1
                                    aplotype = index + 1
                                    i._sample_indexes[sample][0].append(aplotype)
                                    hf_index = variant[3].index(allele)
                                    i._sample_indexes[sample][2].append(variant[7][hf_index])
                                    i._sample_indexes[sample][3].append(variant[8][hf_index])
                                    i._sample_indexes[sample][4].append(variant[9][hf_index])
                                    i._sample_indexes[sample][5].append(variant[5][hf_index])
                        else:
                            i.INFO['AN'] += 1
                            i.samples.append(sample)
                            genotype=[]
                            i._sample_indexes.setdefault(sample,
                                                         [[0], variant[2], variant[7],
                                                          variant[8], variant[9],
                                                          variant[5]])
                            for allele in variant[3]:
                                if allele not in i.ALT:
                                    i.REF.append(variant[1])
                                    i.ALT.append(allele)
                                    i.INFO['AC'].append(1)
                                    i.INFO['AN'] += 1
                                    index = i.ALT.index(allele)
                                    aplotype = index+1
                                    genotype.append(aplotype)
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i.TYPEVAR.append(variant[-1])
                                else:
                                    index = i.ALT.index(allele)
                                    i.INFO['AC'][index] += 1
                                    i.INFO['AN'] += 1
                                    aplotype = index+1
                                    genotype.append(aplotype)
                                    i._sample_indexes[sample][0].append(aplotype)
            else:
                # for homoplasmic variants in a position encountered before
                for i in VCF_RECORDS:
                    if i.POS == variant[0]:
                        for allele in variant[3]:
                            if allele not in i.ALT:
                                i.INFO['AC'].append(1)
                                i.INFO['AN'] += 1
                                i.ALT.append(allele)
                                i.samples.append(sample)
                                index = i.ALT.index(allele)
                                aplotype = index + 1
                                genotype = aplotype
                                i.TYPEVAR.append(variant[-1])
                                if sample in i._sample_indexes:
                                    i._sample_indexes[sample][0].append(genotype)
                                    i._sample_indexes[sample][2].append(variant[7][0])
                                    i._sample_indexes[sample][3].append(variant[8][0])
                                    i._sample_indexes[sample][4].append(variant[9][0])
                                    i._sample_indexes[sample][5].append(variant[5][0])
                                else:
                                    i._sample_indexes.setdefault(sample,
                                                                 [genotype, variant[2],
                                                                  variant[7], variant[8],
                                                                  variant[9], variant[5]])
                                # if a deletion, add a further reference base
                                if type(variant[1]) == type(list()):
                                    i.REF.extend(variant[1])
                                else:
                                    i.REF.append(variant[1])
                            else:
                                index = i.ALT.index(allele)
                                i.INFO['AC'][index] += 1
                                i.INFO['AN'] += 1
                                i.samples.append(sample)
                                aplotype = index + 1
                                genotype = aplotype
                                if sample in i._sample_indexes:
                                    i._sample_indexes[sample][0].append(genotype)
                                    i._sample_indexes[sample][2].append(variant[7][0])
                                    i._sample_indexes[sample][3].append(variant[8][0])
                                    i._sample_indexes[sample][4].append(variant[9][0])
                                    i._sample_indexes[sample][5].append(variant[5][0])
                                else:
                                    i._sample_indexes.setdefault(sample,
                                                                 [genotype, variant[2],
                                                                  variant[7], variant[8],
                                                                  variant[9], variant[5]])

    for r in VCF_RECORDS:
        if len(r.REF) > 1:
            setref = set(r.REF)
            if len(setref) > 1:
                for x in range(len(r.TYPEVAR)):
                    ord = sorted(r.REF, key=lambda x: len(x))
                    if r.TYPEVAR[x] == 'ins':
                        r.ALT[x] = ord[-1] + r.ALT[x][1:]
                    elif r.TYPEVAR[x] == 'del':
                        ndel = len(r.REF[x][1:])
                        altdel = ord[-1][0: (len(ord[-1]) - ndel)]
                        r.ALT[x] = altdel
                    else:
                        r.ALT[x] = r.ALT[x] + ord[-1][1:]
                r.REF = [ord[-1]]
            else:
                r.REF = [r.REF[0]]
    # gets the index of each sample and assign the definitive genotype to each individual
    for index, sample in enumerate(dict_of_dicts.keys()):
        for r in VCF_RECORDS:
            if sample in r.samples:
                r._sample_indexes[sample].append(index)
                if type(r._sample_indexes[sample][0]) == type(list()):
                    genotype = list(map(lambda x: str(x), r._sample_indexes[sample][0]))
                    aplotypes = "/".join(genotype)
                    r._sample_indexes[sample][0] = aplotypes

    # counts also alleles identical to the reference base
    INDEX = OrderedDict()
    for index, sample in enumerate(dict_of_dicts.keys()):
        INDEX.setdefault(sample, index)

    for samples in INDEX:
        for r in VCF_RECORDS:
            if samples not in r._sample_indexes.keys():
                r._sample_indexes.setdefault(samples, [0, INDEX[samples]])
                r.INFO['AN'] += 1
    # VCF header
    header = OrderedDict()
    c = 0
    for x in dict_of_dicts:
        header.setdefault(x, c)
        c += 1

    header = "\t".join(header.keys())

    # writes variant call in the VCF file
    out = BgzfWriter(vcffile, 'w')
    out.write('##fileformat=VCFv4.0\n##reference={}\n'.format(reference))
    out.write('##contig=<ID={},length={}>\n'.format(seq_name, seq_length))
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads covering the REF position">\n')
    out.write('##FORMAT=<ID=HF,Number=3,Type=Float,Description="Heteroplasmy Frequency of variant allele">\n')
    out.write('##FORMAT=<ID=CILOW,Number=3,Type=Float,Description="Value defining the lower limit of the confidence interval of the heteroplasmy fraction">\n')
    out.write('##FORMAT=<ID=CIUP,Number=3,Type=Float,Description="Value defining the upper limit of the confidence interval of the heteroplasmy fraction">\n')
    out.write('##FORMAT=<ID=SDP,Number=2,Type=String,Description="Strand-specific read depth of the ALT allele">\n')
    out.write('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count in genotypes">\n')
    out.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + header + '\n')

    for position in sorted(present_pos):
        for r in VCF_RECORDS:
            if position == r.POS:
                if len(r.INFO['AC'])>1:
                    alleles = list(map(lambda x: str(x), r.INFO['AC']))
                    alleles = ",".join(alleles)
                    AC = 'AC=' + alleles
                    AN = 'AN=' + str(r.INFO['AN'])
                else:
                    AC = 'AC=' + str(r.INFO['AC'][0])
                    AN = 'AN=' + str(r.INFO['AN'])
                samples_per_position = []
                r._sample_indexes = sorted(r._sample_indexes.items(),
                                           key=lambda x: x[1][-1])
                for items in r._sample_indexes:
                    if len(items[1]) > 2:
                        if len(items[1][2]) > 1:
                            het = list(map(lambda x: str(x), items[1][2]))
                            heteroplasmy = ",".join(het)
                            CILOW = list(map(lambda x: str(x), items[1][3]))
                            CIUP = list(map(lambda x: str(x), items[1][4]))
                            confidence_interval_low = ",".join(CILOW)
                            confidence_interval_up = ",".join(CIUP)
                            strandp = []
                            for x in items[1][5]:
                                if isinstance(x, list):
                                    strandp.append(x[0])
                                else:
                                    strandp.append(str(x))
                            SDP = ",".join([str(k) for k in strandp])
                            individual = (str(items[1][0]) + ':' +
                                          str(items[1][1]) + ':' +
                                          str(heteroplasmy) + ':' +
                                          str(confidence_interval_low) + ':' +
                                          str(confidence_interval_up) + ':' +
                                          str(SDP))
                        else:
                            heteroplasmy = str(items[1][2][0])
                            confidence_interval_low = str(items[1][3][0])
                            confidence_interval_up = str(items[1][4][0])
                            if isinstance(items[1][5][0], list):
                                SDP = items[1][5][0][0]
                            else:
                                SDP = str(items[1][5][0])
                            individual = (str(items[1][0]) + ':' +
                                          str(items[1][1]) + ':' +
                                          str(heteroplasmy) + ':' +
                                          str(confidence_interval_low) + ':' +
                                          str(confidence_interval_up) + ':' +
                                          str(SDP))
                        samples_per_position.append(individual)
                    else:
                        individual = str(items[1][0])
                        samples_per_position.append(individual)
                samples = "\t".join(samples_per_position)
                if len(r.ALT) > 1:
                    var = ",".join(r.ALT)
                    out.write(r.CHROM + '\t' + str(r.POS) + '\t' + r.ID + '\t' +
                              r.REF[0] + '\t' + var + '\t' + r.QUAL + '\t' +
                              r.FILTER + '\t' + AC + ';' + AN + '\t' +
                              r.FORMAT + '\t' + samples + '\n')
                else:
                    out.write(r.CHROM + '\t' + str(r.POS) + '\t' + r.ID + '\t' +
                              r.REF[0] + '\t' + r.ALT[0] + '\t' + r.QUAL + '\t' +
                              r.FILTER + '\t' + AC + ';' + AN + '\t' +
                              r.FORMAT + '\t' + samples + '\n')
    out.close()
    return VCF_RECORDS


def FASTAoutput(Consensus, mtDNAseq, names):
    path = os.getcwd()
    fasta_dict2 = {}
    for name2 in names:
        fasta_dict2[name2] = []
    for name2 in fasta_dict2:
        for i in range(len(mtDNAseq)):
            index = i
            val = (index, mtDNAseq[i])
            fasta_dict2[name2].append(val)
    for name2 in Consensus:
        for variants in Consensus[name2]:
            if variants[-1] == 'ins':
                var_pos = (int(variants[0]) - 1) + (float(len(variants[1])) / 10)
                tupla = (var_pos, variants[1])
                fasta_dict2[name2].append(tupla)
            elif variants[-1] == 'del':
                for x in variants[1]:
                    for j in fasta_dict2[name2]:
                        if x == j[0]:
                            fasta_dict2[name2].remove(j)
            else:
                for j in fasta_dict2[name2]:
                    if variants[0] == j[0]:
                        index = fasta_dict2[name2].index(j)
                        fasta_dict2[name2][index] = (variants[0], variants[1])
    for name2 in fasta_dict2:
        for dirname, dirnames, filenames in os.walk('.'):
            for subdirname in dirnames:
                if subdirname.startswith('OUT') and subdirname == names[name2]:
                    fasta_dir = glob.glob(os.path.join(path + '/' + subdirname))[0]
                    fasta_out = open(fasta_dir + '/' + name2 + '.fasta', "w")
                    fasta_out.write('>' + name2 + '_complete_mitochondrial_sequence\n')
                    seq = []
                    fasta_dict2[name2] = sorted(fasta_dict2[name2])
                    for tuples in fasta_dict2[name2]:
                        seq.append(tuples[1][0])
                    seq_def = ''.join(seq)
                    fasta_out.write(seq_def)
                    fasta_out.close()


if __name__ == '__main__':
    print("This script is used only when called by assembleMTgenome.py.")
    pass
