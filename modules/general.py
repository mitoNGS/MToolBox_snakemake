#!/usr/bin/env python

from Bio import SeqIO
import resource, sys, gzip, bz2, re

def memory_usage_resource():
    """
    Function to get memory usage (in MB)
    Source: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    """
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

def s_encoding(s):
    if type(s) == bytes:
        return s.decode("utf-8")
    elif type(s) == str:
        return s

def get_SAM_header(samfile):
    # is the file compressed?
    if samfile.endswith("gz"):
        samhandle = gzip.GzipFile(samfile, mode = 'r')
    elif samfile.endswith("bz2"):
        samhandle = bz2.BZ2File(samfile, mode = 'r')
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

### Functions taken or adapted from assembleMTgenome.py
r=re.compile("#+")
r1=re.compile("""\^.{1}""")
rr=re.compile("[\+\-]{1}[0-9]+")

# variables which should be customisable
mqual=25
clev=0.80
cov=5
glen=10
basename='mtDNAassembly'
sexe='samtools'
sversion=0
crf=0
crc=0
cru=0
pout=0
normb=0
addv=''
addd=''
hf=float(0.8)
tail=5
#

def normS(s,ref):
    c=re.finditer(rr,s)
    sl=list(s)
    cc=[(x.start(),x.end()) for x in c]
    for i in cc:
        n=int(''.join(sl[i[0]+1:i[1]]))
        sl[i[0]:i[1]+n]=['#' for xx in range(len(sl[i[0]:i[1]+n]))]
    ns=''.join(sl)
    ns=ns.replace('#','')
    ss=''
    for i in ns:
        if i in '.,ACGTNacgtN<>*': ss+=i
    return (ss.replace('.',ref)).replace(',',ref)

def nuc(seq):
    d={'A':0,'C':0,'G':0,'T':0,'N':0}
    for i in seq:
        if i in d: d[i]+=1
        else: d['N']+=1
    return d

# dn={'A':'T','T':'A','C':'G','G':'C'}
# def comp(s):
#     ss=''
#     for i in s:
#         if dn.has_key(i): ss+=dn[i]
#         else: ss+='N'
#     return ss

def ff(v,l):
    for i in l:
        x=0
        for j in i:
            if j in v: x+=1
        if x==len(v): return i
    return 0

dIUPAC={'AG':'R','CT':'Y','GC':'S','AT':'W','GT':'K','AC':'M','CGT':'B','AGT':'D','ACT':'H','ACG':'V'}
def getIUPAC(f):
    vv=''.join([i[1] for i in f if i[0]>0])
    k=ff(vv,dIUPAC.keys())
    if k!=0: return dIUPAC[k]
    else: return '#'

def rev(seq):
    d={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    s=''.join([d[x] for x in seq])
    return s[::-1]

def freq(d):
    f=[]
    for i in d:
        try: v=float(d[i])/sum(d.values())
        except: v=0.0
        f.append((v,i))
    f.sort()
    f.reverse()
    maxv=[f[0]]
    for i in f[1:]:
        if i[0]==maxv[0][0]: maxv.append(i)
    if len(maxv)==1:
        if maxv[0][0]>=clev: return maxv[0][1]
        else: return getIUPAC(f)
    elif len(maxv)>1: return getIUPAC(f)

def get_seq_name(fasta):
    mt_genome = SeqIO.index(fasta, 'fasta')
    if len(mt_genome) != 1:
        sys.exit("Sorry, but MToolBox at the moment only accepts single-contig reference mt genomes.")
    for contig, contig_seq in mt_genome.items():
        seq_name = contig
    return seq_name

### FORMAT CONVERTERS

def pileup2mt_table(pileup=None, fasta=None, mt_table=None):
    # generate mt_table backbone
    # get mt sequence data from genome fasta file
    mtdna = {}
    mt_genome = SeqIO.index(fasta, 'fasta')
    if len(mt_genome) != 1:
        sys.exit("Sorry, but MToolBox at the moment only accepts single-contig reference mt genomes.")
    for contig, contig_seq in mt_genome.items():
        for pos, nt in enumerate(contig_seq.seq):
            mtdna[pos+1] = (nt, ['#',(0,0,0,0),0,0.0])

    # open input and output files
    f = open(pileup, 'r')
    mt_table_handle = open(mt_table, 'w')

    # iterate over pileup
    for i in f:
        if i.strip()=='': continue
        l=(i.strip()).split('\t')
        #if l[0]!=mtdna_fasta.split('.')[0]: continue
        pos=int(l[1])
        if len(l) == 6:
            ref,seq,qual=l[2],normS(re.sub(r1,"",l[4]),l[2]),l[5]
            s,q='',0
            for j in range(len(seq)):
                if seq[j] not in '<>*' and ord(qual[j])-33 >= mqual:
                    s+=seq[j].upper()
                    q+=(ord(qual[j])-33)
            try: mq=float(q)/len(s)
            except: mq=0.0
            dnuc=nuc(s)
            mfreq=freq(dnuc)
            lnuc=(dnuc['A'],dnuc['C'],dnuc['G'],dnuc['T'])
            cnuc='#'
            if len(s) >= cov: cnuc=mfreq
            #print pos,cnuc,s,dnuc
            mtdna[pos][1][0]=cnuc
            mtdna[pos][1][1]=lnuc
            mtdna[pos][1][2]=len(s)
            mtdna[pos][1][3]=mq
        else:
            mtdna[pos][1][0]='#'
    f.close()

    # write to mt_table
    aseq = ''
    mt_table_handle.write('Position\tRefNuc\tConsNuc\tCov\tMeanQ\tBaseCount(A,C,G,T)\n')
    assb,totb=0,0
    cop=0
    maxCval=1
    for i in range(len(mtdna)):
        #print i+1, mtdna[i+1]
        line=[str(i+1),mtdna[i+1][0],mtdna[i+1][1][0],str(mtdna[i+1][1][2]),"%.2f" %(mtdna[i+1][1][3]),str(mtdna[i+1][1][1])]
        mt_table_handle.write('\t'.join(line)+'\n')
        #aseq+=mtdna[i+1][1][0]
        # if variant is not #, contigs will have reference, otherwise the # that will be subsequently substituted with N
        if mtdna[i+1][1][0] !='#':
            aseq+=mtdna[i+1][0]
        else:
            aseq+=mtdna[i+1][1][0]
        totb+=1
        if mtdna[i+1][1][0] !='#':
            assb+=1
            cop+=mtdna[i+1][1][2]
            # track.append('chrRSRS %i %i %i\n' %(i,i+1,mtdna[i+1][1][2]))
            if mtdna[i+1][1][2] > maxCval: maxCval=mtdna[i+1][1][2]

    mt_table_handle.close()

### End of functions taken from assembleMTgenome



def sam2fastq(samfile = None, outmt1 = None, outmt2 = None, outmt = None):
    print('Extracting FASTQ from SAM...')
    mtoutsam=samfile
    mtoutfastq=gzip.GzipFile(outmt, 'wb')
    mtoutfastq1=gzip.GzipFile(outmt1, 'wb')
    mtoutfastq2=gzip.GzipFile(outmt2, 'wb')
    f=gzip.GzipFile(mtoutsam, 'rb')
    dics = {}
    c = 0
    for i in f:
        i = i.decode("utf-8")
        c += 1
        if c % 100000 == 0:
            print("{} SAM entries processed.".format(c))
        if i.strip()=='' or i.startswith('@'):
            continue
        l=(i.strip()).split('\t')
        if l[2]=='*': continue
        if len(dics) == 0:
            dics[l[0]]=[l]
        else:
            if l[0] in dics:
                dics[l[0]].append(l)
            else:
                # check if reads go to single or paired end file
                # check if each read in a pair goes to R1 or R2
                if len(dics) != 1:
                    sys.exit("read data not valid: {}".format(dics))
                k = [key for key in dics][0]
                ll=dics[k]
                if len(ll)==1:
                    strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
                    if strand==16: seq,qual=rev(seq),qual[::-1]
                    entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
                    mtoutfastq.write(entry.encode("utf-8"))
                else:
                    strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
                    if strand==16: seq,qual=rev(seq),qual[::-1]
                    entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
                    mtoutfastq1.write(entry.encode("utf-8"))
                    strand,seq,qual=int(ll[1][1]) & 16,ll[1][9],ll[1][10]
                    if strand==16: seq,qual=rev(seq),qual[::-1]
                    entry='\n'.join(['@'+ll[1][0],seq,'+',qual])+'\n'
                    mtoutfastq2.write(entry.encode("utf-8"))
                # create new dics with new read ID
                dics = {l[0] : [l]}
    mtoutfastq.close()
    mtoutfastq1.close()
    mtoutfastq2.close()
