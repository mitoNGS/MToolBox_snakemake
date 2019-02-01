import pandas as pd
import os, re, sys, time, gzip, bz2
from Bio import SeqIO
import resource
import numpy as np
#import sqlite3
from sqlalchemy import create_engine
from modules.mtVariantCaller import *

localrules: bam2pileup, index_genome, pileup2mt_table, make_single_VCF

#shell.prefix("module load gsnap; ")
# fields: sample  ref_genome_mt   ref_genome_n
analysis_tab = pd.read_table("data/analysis.tab", sep = "\t", comment='#')
#print(analysis_tab)
reference_tab = pd.read_table("data/reference_genomes.tab", sep = "\t", comment='#').set_index("ref_genome_mt", drop=False)
#print(reference_tab)

configfile: "config.yaml"
#clusterfile: "cluster."
res_dir = config["results"]
map_dir = config["map_dir"]
log_dir = config["log_dir"]
gmap_db_dir = config["map"]["gmap_db_dir"]
# res_dir = "results"
# map_dir = "map"

def memory_usage_resource():
    """
    Function to get memory usage (in MB)
    Source: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    """
    import resource
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

def read_sam_file_only_readID_chunks_intoSQL(samfile, n_occurrences = 1, chunksize = 100000):
    """
    Read a SAM file, then keep a list of IDs of reads occurring <n_occurrences> times.
    In the specific case of outS.sam and outP.sam files, these are the IDs of the reads we want
    to keep in the OUT.sam file.
    Load the entries in a SQL db
    """
    n = n_occurrences

    # Create in-memory SQLite db
    engine = create_engine('sqlite://', echo=False)
    # samfile = path/to/out.sam --> db_name = out
    db_name = samfile.split('/')[-1].split('.')[0]

    # Read the SAM file in chunks
    header_lines, comment_count = get_SAM_header(samfile)
    # function that reads a samfile and skips rows with unaligned reads
    t = pd.read_table(samfile, \
                      sep = '\t', \
                      skiprows=comment_count, \
                      chunksize=chunksize, \
                      usecols=[0,2], \
                      names = ['readID', 'RNAME'], \
                      compression="infer")

    for chunk in t:
        elapsed = time.time()
        #print("chunk")
        chunk = chunk.query('RNAME != "*"')
        chunk = chunk.drop(columns=['RNAME'])
        chunk.to_sql(db_name, con=engine, if_exists="append")
        print("{} seconds, memory: {} MB".format(time.time()-elapsed, memory_usage_resource()))

    return engine

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

dn={'A':'T','T':'A','C':'G','G':'C'}
def comp(s):
    ss=''
    for i in s:
        if dn.has_key(i): ss+=dn[i]
        else: ss+='N'
    return ss

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

def get_single_vcf_files(df, ref_genome_mt = None):
    ref_genome_mt = ref_genome_mt
    # ref_genome_mt = ref_genome_mt
    # ref_genome_n = ref_genome_n
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
        # "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/vcf.vcf"
            outpaths.append("{}/OUT_{}_{}_{}/vcf.vcf".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n")))
    return outpaths

def get_out_files(df, res_dir="results", map_dir="map"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{}/OUT_{}_{}_{}/{}/OUT.sam".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n"), map_dir))
    return outpaths

def get_vcf_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{}/OUT_{}_{}_{}/vcf.vcf".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n")))
    return outpaths

def get_genome_files(df, ref_genome_mt, field):
    return expand(df.loc[ref_genome_mt, field])

def get_mt_genomes(df):
    return list(set(df['ref_genome_mt']))

def get_mt_fasta(df, ref_genome_mt, field):
    return df.loc[df['ref_genome_mt'] == ref_genome_mt, field][0]

def get_other_fields(df, ref_genome_mt, field):
    return list(set(df.loc[df['ref_genome_mt'] == ref_genome_mt, field]))

def rev(seq):
    d={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    s=''.join([d[x] for x in seq])
    return s[::-1]

def sam2fastq(samfile = None, outmt1 = None, outmt2 = None, outmt = None):
    print('Extracting FASTQ from SAM...')
    #mtoutsam=os.path.join(folder,'samfile')
    mtoutsam=samfile
    mtoutfastq=gzip.GzipFile(outmt, 'w')
    mtoutfastq1=gzip.GzipFile(outmt1, 'w')
    mtoutfastq2=gzip.GzipFile(outmt2, 'w')
    f=gzip.GzipFile(mtoutsam, 'r')
    dics = {}
    c = 0
    for i in f:
        c += 1
        if c % 100000 == 0:
            print("{} SAM entries processed.".format(c))
        # original version
        # if i.strip()=='': continue
        if i.strip()=='' or i.startswith('@'):
            continue
        l=(i.strip()).split('\t')
        #print(l[0])
        if l[2]=='*': continue
        if len(dics) == 0:
            dics[l[0]]=[l]
            #print(dics)
        else:
            if l[0] in dics:
                dics[l[0]].append(l)
                #print(dics)
            else:
                # check if reads go to single or paired end file
                # check if each read in a pair goes to R1 or R2
                if len(dics) != 1:
                    sys.exit("read data not valid: {}".format(dics))
                k = [key for key in dics][0]
                ll=dics[k]
                #print(ll)
                if len(ll)==1:
                    strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
                    if strand==16: seq,qual=rev(seq),qual[::-1]
                    entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
                    mtoutfastq.write(entry)
                    #single.append(entry)
                else:
                    strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
                    if strand==16: seq,qual=rev(seq),qual[::-1]
                    entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
                    mtoutfastq1.write(entry)
                    #pair1.append(entry)
                    strand,seq,qual=int(ll[1][1]) & 16,ll[1][9],ll[1][10]
                    if strand==16: seq,qual=rev(seq),qual[::-1]
                    entry='\n'.join(['@'+ll[1][0],seq,'+',qual])+'\n'
                    mtoutfastq2.write(entry)
                    #pair2.append(entry)
                # create new dics with new read ID
                dics = {l[0] : [l]}
                #print(dics)
    mtoutfastq.close()
    mtoutfastq1.close()
    mtoutfastq2.close()

def read_sam_file_extended(samfile):
    # function that reads a samfile and skips rows with unaligned reads
    t = pd.read_table(samfile, sep = '\t', comment='@', usecols=[0,2,3], names = ['readID', 'RNAME', 'POS'])
    t = t.loc[t['RNAME'] != '*']
    return t

def read_sam_file(samfile):
    # function that reads a samfile and skips rows with unaligned reads
    t = pd.read_table(samfile, sep = '\t', comment='@', usecols=[0,2], names = ['readID', 'RNAME'])
    t = t.loc[t['RNAME'] != '*']
    return t

def filter_alignments(outmt = None, outS = None, outP = None, OUT = None, ref_mt_fasta = None):
    print("Processing {}".format(outS))
    outS_sql = read_sam_file_only_readID_chunks_intoSQL(outS)
    print("Processing {}".format(outP))
    outP_sql = read_sam_file_only_readID_chunks_intoSQL(outP)

    #print(outS_sql.table_names())
    #print(outP_sql.table_names())
    #outS_sql.execute("SELECT readID from outS GROUP BY readID")
    #pd.read_sql_query("SELECT readID, COUNT(readID) FROM outS HAVING COUNT(readID) == 1 GROUP BY readID", outS_sql)
    good_reads = pd.concat([pd.read_sql_query("SELECT readID FROM outS GROUP BY readID HAVING COUNT(*) == 1", outS_sql), \
                            pd.read_sql_query("SELECT readID FROM outP GROUP BY readID HAVING COUNT(*) == 2", outP_sql)])
    print("Total reads to extract alignments of: {}".format(len(good_reads)))

    samfile = outmt
    tc = pd.read_table(samfile, \
                       sep = '\t', \
                       skiprows=get_SAM_header(samfile)[1], \
                       chunksize=100000, \
                       header=None, \
                       engine="python", \
                       names=["readID", "FLAG", "RNAME"] + list("QWERTYUIOPASDFGHJK"), \
                       compression="infer")

    # open OUT.sam file and write SAM header from outS.sam (outP would be the same).
    OUT_uncompressed = OUT.replace(".gz", "")
    f = open(OUT_uncompressed, 'w')
    sss = gzip.open(outS, 'rb')
    l = sss.readline().decode("utf-8")
    while l[0] == "@":
        if l.startswith("@PG") == False:
            f.write(l)
        l = sss.readline().decode("utf-8")

    f.close()

    n_extracted_alignments = 0
    #D = pd.DataFrame()
    for chunk in tc:
        chunk = chunk.query('RNAME != "*"')
        OUT_chunk = pd.merge(chunk, good_reads, how="inner", on="readID")
        #c_se = pd.merge(chunk, se, how="inner", on="readID")
        #c_pe = pd.merge(chunk, pe, how="inner", on="readID")
        #print(len(c))
        n_extracted_alignments += len(OUT_chunk)
        # Append alignments to OUT.sam
        OUT_chunk.to_csv(OUT_uncompressed, mode="a", header=False, sep="\t", index=False)
        #n_extracted_reads += len(c_se)
        #n_extracted_reads += len(c_pe)
        #D = pd.concat([D, OUT_chunk])
    
    print("Compressing OUT.sam file")
    os.system("gzip {}".format(OUT_uncompressed))
    print("Total alignments extracted: {}".format(n_extracted_alignments))

# def filter_alignments_old(outmt = None, outS = None, outP = None, OUT = None, ref_mt_fasta = None):
#     print("Memory usage at the beginning of the function: {} MB".format(memory_usage_resource()))
#     ref_mt_fasta = SeqIO.parse(ref_mt_fasta, 'fasta')
#     sig=1
#     pai=1
#     print('Reading Results...')
#     if sig:
#         hgoutsam = outS
#         dicsingle={}
#         f = open(hgoutsam, 'r')
#         x = 0
#         for i in f:
#             if i.strip()=='': continue
#             l=(i.strip()).split('\t')
#             if l[2]=='*': continue # the read is not mapped
#             # keeping multiple mappings
#             if l[0] in dicsingle:
#                 dicsingle[l[0]].append(l)
#             else:
#                 dicsingle[l[0]]=[l]
#             x += 1
#             if x%100000 == 0:
#                 print("Memory usage after {} reads: {} MB".format(x, memory_usage_resource()))
#         f.close()
#     if pai:
#         hgoutsam2 = outP
#         dicpair={}
#         f=open(hgoutsam2)
#         x = 0
#         for i in f:
#             if i.strip()=='': continue
#             l=(i.strip()).split('\t')
#             if l[2]=='*': continue
#             if l[0] in dicpair:
#                 dicpair[l[0]].append(l)
#             else:
#                 dicpair[l[0]]=[l]
#             x += 1
#             if x%100000 == 0:
#                 print("Memory usage after {} reads: {} MB".format(x, memory_usage_resource()))
#         f.close()
# 
#     #print('Extracting FASTQ from SAM...')
#     mtoutsam = outmt
#     dics={}
#     f=open(mtoutsam)
#     for i in f:
#         if i.strip()=='' or i.startswith('@'): continue
#         l=(i.strip()).split('\t')
#         if l[2]=='*': continue
#         if l[0] in dics: dics[l[0]].append(l)
#         else: dics[l[0]]=[l]
#     f.close()
# 
#     finalsam = OUT
#     out=open(finalsam,'w')
#     for entry in ref_mt_fasta:
#         out.write("@SQ    SN:{}	LN:{}\n".format(entry.id, len(entry)))
#     #out.write("@SQ	SN:%s	LN:16569\n" % gsnap_db)
#     out.write("@RG	ID:sample	PL:sample	PU:sample	LB:sample	SM:sample\n")
# 
#     print('Filtering reads...')
#     for i in dics:
#         ll=dics[i]
#         if len(ll)==1: # if the read has one mapping I assume it's SE
#             if i in dicsingle:
#                 r=dicsingle[i] # i is a list of lists (splitted sam lines)
#                 if len(r)==1:
#                     # check if read aligned on MT when aligned against nuclear+MT
#                     # fields checked: RNAME, POS
#                     if r[0][2]==ll[0][2] and ll[0][3]==r[0][3]:
#                         #good.append('\t'.join(ll[0])+'\n')
#                         out.write('\t'.join(ll[0])+'\n')
#             else:
#                 out.write('\t'.join(ll[0])+'\n')
#         else:
#             if i in dicpair:
#                 r=dicpair[i]
#                 if len(r) == 2:
#                     if r[0][2]==ll[0][2] and ll[0][3]==r[0][3] and r[1][2]==ll[1][2] and ll[1][3]==r[1][3]:
#                         out.write('\t'.join(ll[0])+'\n')
#                         out.write('\t'.join(ll[1])+'\n')
#             else:
#                 out.write('\t'.join(ll[0])+'\n')
#                 out.write('\t'.join(ll[1])+'\n')
#     out.close()
# 
#     print('Outfile saved on %s.' %(finalsam))
#     print('Done.')

def gsnap_inputs(sample, read_type):
    # https://stackoverflow.com/questions/6930982/how-to-use-a-variable-inside-a-regular-expression
    ### This regex matches typical names in illumina sequencing, eg 95191_TGACCA_L002_R2_001.fastq.gz
    #read_file_regex = re.escape(sample) + r'_[\D]{6}_L001_R' + read_type + r'_001.fastq.gz'
    ### This is for more general cases, seems to work
    read_file_regex = re.escape(sample) + r'[_.][\S]*R' + read_type + r'[\S]*fastq.gz'
    read_file = [f for f in os.listdir('./data/reads/') if re.match(read_file_regex, f)]
    if len(read_file) > 1:
        sys.exit("Ambiguous name in read files.")
    elif len(read_file) == 0:
        sys.exit("No read files found for sample: {}".format(sample))
    return 'data/reads/' + read_file[0]

wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in list(set(analysis_tab['sample']))]),
    ref_genome_mt = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_mt']))]),
    ref_genome_n = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_n']))]),

outpaths = get_mt_genomes(analysis_tab)

target_inputs = [
    outpaths ]

# rule all:
#     input:
#         mt_vcf = expand("results/vcf/{ref_genome_mt}.vcf", ref_genome_mt = get_mt_genomes(analysis_tab)),

rule all:
    input:
        get_vcf_files(analysis_tab)
        # single_vcf = lambda wildcards: expand("results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/vcf.vcf",
        #                                         sample = get_other_fields(analysis_tab, wildcards.ref_genome_mt, "sample"),
        #                                         ref_genome_mt = get_mt_genomes(analysis_tab),
        #                                         #ref_genome_mt = wildcards.ref_genome_mt,
        #                                         ref_genome_n = get_other_fields(analysis_tab, wildcards.ref_genome_mt, "ref_genome_n"))

rule make_mt_gmap_db:
    input:
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}", \
                            ref_genome_mt_file = get_genome_files(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
    output:
        gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].split(".")[0]
    message: "Generating gmap db for mt genome: {input.mt_genome_fasta}"
    shell:
        """
        module load gsnap
        gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -s numeric-alpha {input.mt_genome_fasta}
        """

rule make_mt_n_gmap_db:
    input:
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}", \
                            ref_genome_mt_file = get_genome_files(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file")),
        n_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_n_file}", \
                            ref_genome_n_file = get_genome_files(reference_tab, wildcards.ref_genome_mt, "ref_genome_n_file"))
    output:
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome",
		mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    params:
        #mt_n_fasta = lambda wildcards: "data/genomes/{}_{}.fasta".format(wildcards.ref_genome_mt, wildcards.ref_genome_n),
        #mt_n_fasta = lambda wildcards, input: "{}_{}.fasta".format(wildcards.ref_genome_mt, os.path.split(input.n_genome_fasta)[1].split(".")[0]),
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].split(".")[0]
    message: "Generating gmap db for mt + n genome: {input.mt_genome_fasta},{input.n_genome_fasta}"
    shell:
        """
        module load gsnap
        cat {input.mt_genome_fasta} {input.n_genome_fasta} > {params.mt_n_fasta}
        gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -s numeric-alpha {params.mt_n_fasta}
        # rm {input.mt_genome_fasta}_{input.n_genome_fasta}.fasta
        """

seq_type = "pe"

rule map_MT_PE_SE:
    input:
        R1 = lambda wildcards: gsnap_inputs("{sample}".format(sample=wildcards.sample), "1"),
        R2 = lambda wildcards: gsnap_inputs("{sample}".format(sample=wildcards.sample), "2"),
        gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    output:
        outmt_sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        gmap_db = lambda wildcards: wildcards.ref_genome_mt,
        RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample'
    log:
        log_dir + "/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/logmt.txt"
    threads:
        config["map"]["gmap_threads"]
    message: "Mapping reads for sample {wildcards.sample} to {wildcards.ref_genome_mt} mt genome"
    run:
        if seq_type == "pe":
            print("PE mode")
            shell("module load gsnap")
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} | gzip -c - > {output.outmt_sam} 2> {log}")
        if seq_type == "se":
            print("SE mode")
            shell("module load gsnap")
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} | gzip -c - > {output.outmt_sam} 2> {log}")
        elif seq_type == "both":
            print("PE + SE mode")
            shell("module load gsnap")
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} {input[2]} | gzip -c - > {output.outmt_sam} 2> {log}")

rule sam2fastq:
    input:
        outmt_sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt.sam.gz"
    output:
        outmt1 = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt1.fastq.gz",
        outmt2 = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt2.fastq.gz",
        outmt = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt.fastq.gz",
        #log = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/sam2fastq.done"
    # version:
    #     subprocess.getoutput(
    #         "picard SamToFastq --version"
    #         )
    message:
        "Converting SAM files to FASTQ"
    run:
        sam2fastq(samfile = input.outmt_sam, outmt1 = output.outmt1, outmt2 = output.outmt2, outmt = output.outmt)
    # shell:
    #     """
    #     picard SamToFastq \
    #         I={input.outmt_sam} \
    #         FASTQ={output.outmt1} \
    #         SECOND_END_FASTQ={output.outmt2} \
    #         UNPAIRED_FASTQ={output.outmt}
    #     """

rule map_nuclear_MT_SE:
    input:
        outmt = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt.fastq.gz",
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome"
    output:
        outS = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outS.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].split(".")[0]
        #gsnap_db = config['map']['gsnap_n_mt_db']
    threads:
        config["map"]["gmap_remap_threads"]
    # version:
    #     subprocess.getoutput(
    #       gsnap --version
    #       )
    #     """
    log:
        logS = log_dir + "/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/logS.sam"
    message:
        "Mapping onto complete human genome (nuclear + mt)... SE reads"
    run:
        shell("module load gsnap")
        #shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt} > {output.outS} 2> {log.logS}")
        shell("if [[ -s {input.outmt} ]]; then gsnap -D {params.gmap_db_dir} -d {params.gmap_db} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt} | gzip -c - > {output.outS} 2> {log.logS}; else touch {output.outS}; fi")
    # shell:
    #     """
    #     #module load gsnap
    #     if [[ -s {input.outmt} ]]
    #     then
    #         gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -A sam --nofails --query-unk-mismatch=1 -O -t {threads} -o {output.outS} {input.outmt} &> {log.logS}
    #     else
    #         touch {output.outS}
    #     fi
    #     """

rule map_nuclear_MT_PE:
    input:
        outmt1 = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt1.fastq.gz",
        outmt2 = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt2.fastq.gz",
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome"
    output:
        outP = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outP.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].split(".")[0]
        #gsnap_db = config['map']['gsnap_n_mt_db']
    threads:
        config["map"]["gmap_remap_threads"]
    # version:
    #     subprocess.getoutput(
    #       gsnap --version
    #       )
    #     """
    log:
        logP = log_dir + "/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/logP.sam"
    message:
        "Mapping onto complete human genome (nuclear + mt)... PE reads"
    shell:
        """
        #module load gsnap
        if [[ -s {input.outmt1} ]]
        then
            gsnap -D {params.gmap_db_dir} -d {params.gmap_db} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt1} {input.outmt2} | gzip -c - > {output.outP} 2> {log.logP}
        else
            touch {output}
        fi
        """

rule filtering_mt_alignments:
    input:
        outmt = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outmt.sam.gz",
        outS = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outS.sam.gz",
        outP = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/outP.sam.gz"
    output:
        sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.sam.gz"
        #sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.sam"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
        #outdir = lambda wildcards, output: os.path.split(output.sam)[0]
    threads: 1
    message: "Filtering alignments in file {input.outmt} by checking alignments in {input.outS} and {input.outP}"
    run:
        filter_alignments(outmt = input.outmt, \
                          outS = input.outS, \
                          outP = input.outP, \
                          OUT = output.sam, \
                          ref_mt_fasta = params.ref_mt_fasta)

rule sam2bam:
    input:
        sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.sam.gz",
        #n_genome_fasta = "data/genomes/{ref_genome_n_file}"
    output:
        "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.bam"
    message: "Converting {input.sam} to {output}"
    shell:
        """
        zcat {input.sam} | samtools view -b -o {output} -
        """

rule sort_bam:
    input:
        bam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.bam"
    output:
        sorted_bam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT-sorted.bam"
    message: "Sorting {input.bam} to {output.sorted_bam}"
    shell:
        """
        samtools sort -o {output.sorted_bam} -T ${{TMP}} {input.bam}
        """

rule index_genome:
    input:
        mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
    message: "Indexing {input.genome_fasta} with samtools faidx"
    shell:
        """
        samtools faidx {input.mt_n_fasta}
        """

rule bam2pileup:
    input:
        sorted_bam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT-sorted.bam",
        genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
    output:
        pileup = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/OUT-sorted.pileup"
    params:
        genome_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    message: "Generating pileup {output.pileup} from {input.sorted_bam}"
    shell:
        """
        samtools mpileup -B -f {params.genome_fasta} -o {output.pileup} {input.sorted_bam}
        """

rule pileup2mt_table:
    input:
        pileup = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/OUT-sorted.pileup"
    output:
        mt_table = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/OUT-mt_table.txt"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
    message: "Generating mt_table {output.mt_table} from {input.pileup}, ref mt: {params.ref_mt_fasta}"
    run:
        pileup2mt_table(pileup=input.pileup, fasta=params.ref_mt_fasta, mt_table=output.mt_table)

rule make_single_VCF:
    input:
        sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.sam.gz",
        mt_table = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/OUT-mt_table.txt"
    output:
        single_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/vcf.vcf"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
    message: "Processing {input.sam} to get VCF {output.single_vcf}"
    run:
        # function (and related ones) from mtVariantCaller
        vcf_dict = mtvcf_main_analysis(sam_file = input.sam, mtable_file = input.mt_table, name2 = wildcards.sample)
        # ref_genome_mt will be used in the VCF descriptive field
        # seq_name in the VCF data
        seq_name = get_seq_name(params.ref_mt_fasta)
        VCFoutput(vcf_dict, reference = wildcards.ref_genome_mt, seq_name = seq_name, vcffile = output.single_vcf)

# rule make_VCF:
#     input:
#         single_vcf = lambda wildcards: expand("results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/vcf.vcf",
#                                                 sample = get_other_fields(analysis_tab, wildcards.ref_genome_mt, "sample"),
#                                                 ref_genome_mt = wildcards.ref_genome_mt,
#                                                 ref_genome_n = get_other_fields(analysis_tab, wildcards.ref_genome_mt, "ref_genome_n"))
#     output:
#         "results/vcf/{ref_genome_mt}.vcf"
#     message: "Merging VCFs: {input.single_vcf} into file {output}"
#     shell:
#         """
#         # do something
#         cmd {input} {output}
#         """
