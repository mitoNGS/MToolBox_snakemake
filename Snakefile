import pandas as pd
import os, re, sys, time, gzip, bz2, subprocess
from Bio import SeqIO
import resource
import numpy as np
#import sqlite3
from sqlalchemy import create_engine
from modules.mtVariantCaller import *
from modules.BEDoutput import *

#localrules: bam2pileup, index_genome, pileup2mt_table, make_single_VCF
localrules: index_genome, merge_VCF, index_VCF, dict_genome

#shell.prefix("module load gsnap; ")
# fields: sample  ref_genome_mt   ref_genome_n
analysis_tab = pd.read_table("data/analysis.tab", sep = "\t", comment='#')
#print(analysis_tab)
reference_tab = pd.read_table("data/reference_genomes.tab", sep = "\t", comment='#').set_index("ref_genome_mt", drop=False)
datasets_tab = pd.read_table("data/datasets.tab", sep = "\t", comment='#')

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

def read_sam_file_only_readID_chunks_intoSQL(samfile, n_occurrences = 1, chunksize = 100000, table_name = "outS", ext = ".sam.gz"):
    """
    Read a SAM file, then keep a list of IDs of reads occurring <n_occurrences> times.
    In the specific case of outS.sam and outP.sam files, these are the IDs of the reads we want
    to keep in the OUT.sam file.
    Load the entries in a SQL db
    """
    n = n_occurrences

    # Create in-memory SQLite db
    engine = create_engine('sqlite://', echo=False)
    # samfile = path/to/out.sam --> table_name = out
    table_name = table_name
    #table_name = samfile.split('/')[-1].replace(ext, "").upper().replace("-", "_")

    # Read the SAM file in chunks
    header_lines, comment_count = get_SAM_header(samfile)
    # function that reads a samfile and skips rows with unaligned reads
    t = pd.read_table(samfile, \
                      sep = '\t', \
                      skiprows=comment_count, \
                      chunksize=chunksize, \
                      usecols=[0,2], \
                      names = ['readID', 'RNAME'], \
                      compression="infer", index_col = False)

    for chunk in t:
        elapsed = time.time()
        #print("chunk")
        chunk = chunk.query('RNAME != "*"')
        chunk = chunk.drop(columns=['RNAME'])
        chunk.to_sql(table_name, con=engine, if_exists="append")
        print("{} seconds, memory: {} MB".format(time.time()-elapsed, memory_usage_resource()))

    return engine, table_name

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
        outpaths.append("{results}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf".format(results = res_dir, \
                                                                                                                    sample = getattr(row, "sample"), \
                                                                                                                    ref_genome_mt = getattr(row, "ref_genome_mt"), \
                                                                                                                    ref_genome_n = getattr(row, "ref_genome_n")))
    return outpaths

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

#sorted_bams = lambda.wildcards: get_sample_bamfiles(df, res_dir="results", sample = wildcards.sample, ref_genome_mt = wildcards.ref_genome_mt, ref_genome_n = wildcards.ref_genome_n)
#"results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT.bam"
def get_sample_bamfiles(df, res_dir="results", sample = None, ref_genome_mt = None, ref_genome_n = None):
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "sample") == sample:
            bam_file = getattr(row, "R1").replace("_R1_001.fastq.gz", "_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam".format(ref_genome_mt = ref_genome_mt, ref_genome_n = ref_genome_n))
            out_folder = "OUT_{base}".format(base = bam_file.replace("_OUT-sorted.bam", ""))
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

def get_other_fields(df, ref_genome_mt, field):
    return list(set(df.loc[df['ref_genome_mt'] == ref_genome_mt, field]))

def rev(seq):
    d={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    s=''.join([d[x] for x in seq])
    return s[::-1]

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

def filter_alignments(outmt = None, outS = None, outP = None, OUT = None, ref_mt_fasta = None):
    print("Processing {}".format(outS))
    outS_sql, table_name_S = read_sam_file_only_readID_chunks_intoSQL(outS, table_name = "outS")
    print("Processing {}".format(outP))
    outP_sql, table_name_P = read_sam_file_only_readID_chunks_intoSQL(outP, table_name = "outP")

    good_reads_S = pd.read_sql_query("SELECT readID FROM {table_name} GROUP BY readID HAVING COUNT(*) == 1".format(table_name=table_name_S), outS_sql, chunksize = 100000)
    print("SQL query on outS, memory: {} MB".format(memory_usage_resource()))
    good_reads_P = pd.read_sql_query("SELECT readID FROM {table_name} GROUP BY readID HAVING COUNT(*) == 2".format(table_name=table_name_P), outP_sql, chunksize = 100000)
    print("SQL query on outP, memory: {} MB".format(memory_usage_resource()))
    good_reads = pd.DataFrame()
    for c in good_reads_P:
        good_reads = good_reads.append(c)
    print("good_reads_P append, memory: {} MB".format(memory_usage_resource()))
    for c in good_reads_S:
        good_reads = good_reads.append(c)
    print("good_reads_S append, memory: {} MB".format(memory_usage_resource()))
    print("Total reads to extract alignments of: {}".format(len(good_reads)))

    samfile = outmt
    tc = pd.read_table(samfile, \
                       sep = '\t', \
                       skiprows=get_SAM_header(samfile)[1], \
                       chunksize=100000, \
                       header=None, \
                       engine="python", \
                       names=["readID", "FLAG", "RNAME"] + list("QWERTYUIOPASDFGHJK"), \
                       compression="infer", index_col = False)

    # open OUT.sam file and write SAM header from outS.sam (outP would be the same).
    OUT_uncompressed = OUT.replace(".gz", "")
    f = open(OUT_uncompressed, 'w')
    sss = gzip.open(outS, 'rb')
    l = sss.readline().decode("utf-8")
    while l[0] == "@":
        if l.startswith("@PG") == False:
            f.write(l)
        l = sss.readline().decode("utf-8")
    f.write("\t".join(["@RG", "ID:sample", "PL:illumina", "SM:sample"])+"\n")

    f.close()

    n_extracted_alignments = 0
    for chunk in tc:
        chunk = chunk.query('RNAME != "*"')
        OUT_chunk = pd.merge(chunk, good_reads, how="inner", on="readID")
        print("Chunk, memory: {} MB".format(memory_usage_resource()))
        n_extracted_alignments += len(OUT_chunk)
        # Append alignments to OUT.sam
        OUT_chunk.to_csv(OUT_uncompressed, mode="a", header=False, sep="\t", index=False)

    print("Compressing OUT.sam file")
    os.system("gzip {}".format(OUT_uncompressed))
    print("OUT.sam compressed, memory: {} MB".format(memory_usage_resource()))
    print("Total alignments extracted: {}".format(n_extracted_alignments))

def read_datasets_inputs(sample = None, read_type = "1", input_folder="data/reads"):
    # https://stackoverflow.com/questions/6930982/how-to-use-a-variable-inside-a-regular-expression
    ### This regex matches typical names in illumina sequencing, eg 95191_TGACCA_L002_R2_001.fastq.gz
    #read_file_regex = re.escape(sample) + r'_[\D]{6}_L001_R' + read_type + r'_001.fastq.gz'
    ### This is for more general cases, seems to work
    read_file_regex = re.escape(sample) + r'[_.][\S]*R' + read_type + r'[\S]*fastq.gz'
    #print(os.listdir(input_folder))
    read_files = [f for f in os.listdir(input_folder) if re.match(read_file_regex, f)]
    #print(read_files)
    if len(read_files) > 1:
        sys.exit("Ambiguous name in read files.")
    elif len(read_files) == 0:
        sys.exit("No read files found for sample: {}".format(sample))
    return [os.path.join(input_folder, r) for r in read_files]

def fastqc_raw_outputs(datasets_tab, analysis_tab = None, infolder="data/reads", outfolder="results/fastqc_raw", ext=".fastq.gz"):
    fastqc_out = []
    for i,l in datasets_tab.iterrows():
        if l["sample"] in list(analysis_tab["sample"]):
            fastqc_out.append(os.path.join(outfolder, l["R1"].replace(ext, "_fastqc.html")))
            fastqc_out.append(os.path.join(outfolder, l["R2"].replace(ext, "_fastqc.html")))
    return fastqc_out

# results/fastqc_filtered/{sample}_{adapter}_{lane}_R1_fastqc.html
# html_report_R1 = "results/fastqc_filtered/{sample}_{adapter}_{lane}_qc_R1_fastqc.html",
# html_report_R2 = "results/fastqc_filtered/{sample}_{adapter}_{lane}_qc_R2_fastqc.html",
# html_report_U = "results/fastqc_filtered/{sample}_{adapter}_{lane}_qc_U_fastqc.html",
def fastqc_filtered_outputs(datasets_tab, analysis_tab = None, infolder="data/reads", outfolder="results/fastqc_filtered", ext="_001.fastq.gz"):
    fastqc_out = []
    for i,l in datasets_tab.iterrows():
        if l["sample"] in list(analysis_tab["sample"]):
            fastqc_out.append(os.path.join(outfolder, l["R1"].replace("_R1_001.fastq.gz", "_qc_R1_fastqc.html")))
            fastqc_out.append(os.path.join(outfolder, l["R2"].replace("_R2_001.fastq.gz", "_qc_R2_fastqc.html")))
            fastqc_out.append(os.path.join(outfolder, l["R1"].replace("_R1_001.fastq.gz", "_qc_U_fastqc.html")))
    return fastqc_out

# def fastqc_raw_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_raw", ext=".fastq.gz"):
#     # analysis_tab is a pandas df with a column "sample"
#     fastqc_out = []
#     for s in analysis_tab["sample"]:
#         # fastqc_html_1 and fastqc_html_2 could be strings, but
#         # keep them as lists in case one sample has multiple datasets
#         fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
#         fastqc_out.extend(fastq_files_1)
#         fastqc_html_2 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "2", input_folder = infolder)]
#         fastqc_out.extend(fastq_files_2)
#     return fastqc_out
#
# def fastqc_filtered_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_filtered", ext=".fastq.gz"):
#     fastqc_out = []
#     for s in analysis_tab["sample"]:
#         # fastqc_html_1 and fastqc_html_2 could be strings, but
#         # keep them as lists in case one sample has multiple datasets
#         #L = set(os.path.join(infolder))
#         #glob.glob("{in}/{}*R1{}.format()")
#         # expand("results/fastqc_filtered/{sample}")
#         # fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
#         # fastqc_out.extend(fastq_html_1)
#         #######
#         for read_type in ["R1", "R2", "U"]:
#             fastqc_out.append(os.path.join(outfolder, "{sample}.{read_type}_fastqc.html".format(sample=s, read_type=read_type)))
#     return fastqc_out

# def fastqc_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_raw", ext=".fastq.gz", read_types = ["1", "2"]):
#     # keyword default values are for raw reads
#     fastqc_out = []
#     for s in analysis_tab["sample"]:
#         # fastqc_html_1 and fastqc_html_2 could be strings, but
#         # keep them as lists in case one sample has multiple datasets
#         for read_type in read_types:
#             for i in read_datasets_inputs(sample = s, read_type = read_type, input_folder = infolder):
#                 print(i)
#                 print(os.path.split(i)[1])
#                 fastqc_html = [os.path.join(outfolder, os.path.split(i)[1].replace(ext, "_fastqc.html"))]
#             fastqc_out.extend(fastqc_html)
#         # fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
#         # fastqc_out.extend(fastq_html_1)
#         # fastqc_html_2 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "2", input_folder = infolder)]
#         # fastqc_out.extend(fastq_html_2)
#         # fastqc_html_U = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "U", input_folder = infolder)]
#         # fastqc_out.extend(fastq_html_U)
#     print("fastqc_outputs: {}".format(fastqc_out))
#     return fastqc_out


wildcard_constraints:
    sample = '|'.join([re.escape(x) for x in list(set(analysis_tab['sample']))]),
    ref_genome_mt = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_mt']))]),
    ref_genome_n = '|'.join([re.escape(x) for x in list(set(analysis_tab['ref_genome_n']))])

outpaths = get_mt_genomes(analysis_tab)

target_inputs = [
    outpaths ]

# rule all:
#     input:
#         get_genome_vcf_files(analysis_tab),
#         get_bed_files(analysis_tab),
#         fastqc_raw_outputs = expand("results/fastqc_raw/{sample}_{read_type}_fastqc.html", sample=analysis_tab["sample"], read_type = ["R1", "R2"]),
#         fastqc_filtered_outputs = fastqc_filtered_outputs(analysis_tab = analysis_tab, infolder = "data/reads", outfolder = "results/fastqc_filtered", ext = ".fastq.gz"),

rule all:
    input:
        fastqc_raw_outputs(datasets_tab, analysis_tab = analysis_tab),
        fastqc_filtered_outputs(datasets_tab, analysis_tab = analysis_tab),
        get_genome_vcf_files(analysis_tab),
        get_bed_files(analysis_tab),

rule fastqc_raw:
    input:
        R1 = "data/reads/{sample}_{lane}_R1_001.fastq.gz",
        R2 = "data/reads/{sample}_{lane}_R2_001.fastq.gz"
        # R1 = lambda wildcards: read_datasets_inputs(sample="{sample}".format(sample=wildcards.sample), read_type="1", input_folder="data/reads"),
        # R2 = lambda wildcards: read_datasets_inputs(sample="{sample}".format(sample=wildcards.sample), read_type="2", input_folder="data/reads"),
    output:
        html_report_R1 = "results/fastqc_raw/{sample}_{lane}_R1_001_fastqc.html",
        html_report_R2 = "results/fastqc_raw/{sample}_{lane}_R2_001_fastqc.html",
        #html_report_R1 = "results/fastqc_raw/{sample}_R1_fastqc.html",
        #html_report_R2 = "results/fastqc_raw/{sample}_R2_fastqc.html",
        # html_report_R1 = lambda wildcards: "results/fastqc_raw/{outR1}".format(outR1 = os.path.split({input.R1})[1].replace(".fastq.gz", "_fastqc.html")),
        # html_report_R2 = lambda wildcards: "results/fastqc_raw/{outR2}".format(outR2 = os.path.split({input.R2})[1].replace(".fastq.gz", "_fastqc.html"))
        #logFile = os.path.join(config["proj_dirs"]["logs"], "fastqc_raw.log")
    params:
        outDir = "results/fastqc_raw/",
    threads:
        2
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of raw read files {input} with {version}, {wildcards}"
    log:
        "logs/fastqc_raw/{sample}_{lane}.log"
    #conda: "envs/environment.yaml"
    shell:
        """
        mkdir -p {params.outDir}
        fastqc -t {threads} -o {params.outDir} {input} &> {log}
        """

rule make_mt_gmap_db:
    input:
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}", \
                            ref_genome_mt_file = get_genome_files(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
    output:
        gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        # gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].split(".")[0]
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
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}", \
                            ref_genome_mt_file = get_genome_files(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file")),
        n_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_n_file}", \
                            ref_genome_n_file = get_genome_files(reference_tab, wildcards.ref_genome_mt, "ref_genome_n_file"))
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
        out1P = "data/reads_filtered/{sample}_{lane}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{lane}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{lane}_qc_U.fastq.gz",
        # R1 = "data/reads_filtered/{sample}.R1.fastq.gz",
        # R2 = "data/reads_filtered/{sample}.R2.fastq.gz",
        # U = "data/reads_filtered/{sample}.U.fastq.gz"
    output:
        html_report_R1 = "results/fastqc_filtered/{sample}_{lane}_qc_R1_fastqc.html",
        html_report_R2 = "results/fastqc_filtered/{sample}_{lane}_qc_R2_fastqc.html",
        html_report_U = "results/fastqc_filtered/{sample}_{lane}_qc_U_fastqc.html",
        # html_report_R1 = "results/fastqc_filtered/{sample}.R1_fastqc.html",
        # html_report_R2 = "results/fastqc_filtered/{sample}.R2_fastqc.html",
        # html_report_U = "results/fastqc_filtered/{sample}.U_fastqc.html",
    params:
        outDir = "results/fastqc_filtered/"
    threads:
        3
    # version:
    #     subprocess.check_output("fastqc -V", shell=True)
    # message:
    #     "QC of filtered read files {input} with {version}"
    log:
        "logs/fastqc_filtered/{sample}_{lane}.log"
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
        jar_file = config['read_processing']['trimmomatic']['jar_file'],
        mem = config['read_processing']['trimmomatic']['java_vm_mem'],
        options = config['read_processing']['trimmomatic']['options'],
        processing_options = config['read_processing']['trimmomatic']['processing_options'],
        out1P = "data/reads_filtered/{sample}_{lane}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{lane}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{lane}_qc_1U.fastq.gz",
        out2U = "data/reads_filtered/{sample}_{lane}_qc_2U.fastq.gz"
    input:
        R1 = "data/reads/{sample}_{lane}_R1_001.fastq.gz",
        R2 = "data/reads/{sample}_{lane}_R2_001.fastq.gz"
        # R1 = lambda wildcards: read_datasets_inputs(sample = "{sample}".format(sample=wildcards.sample), read_type = "1", input_folder = "data/reads"),
        # R2 = lambda wildcards: read_datasets_inputs(sample = "{sample}".format(sample=wildcards.sample), read_type = "2", input_folder = "data/reads"),
    output:
        out1P = "data/reads_filtered/{sample}_{lane}_qc_R1.fastq.gz",
        out2P = "data/reads_filtered/{sample}_{lane}_qc_R2.fastq.gz",
        out1U = "data/reads_filtered/{sample}_{lane}_qc_U.fastq.gz",
        # out1P = config['proj_dirs']['reads_filtered'] + "/{sample}.R1.fastq.gz",
        # out2P = config['proj_dirs']['reads_filtered'] + "/{sample}.R2.fastq.gz",
        # out1U = config['proj_dirs']['reads_filtered'] + "/{sample}.fastq.gz",
    threads:
        config['read_processing']['trimmomatic']['threads']
    # version:
    #     subprocess.check_output("trimmomatic -version", shell=True)
    message:
        "Filtering read dataset {wildcards.sample}_{wildcards.lane} with Trimmomatic. {wildcards}" # v{version}"
    log:
        log_dir + "/trimmomatic/{sample}_{lane}_trimmomatic.log"
    #conda: "envs/environment.yaml"
    shell:
        """
        {params.java_cmd} -Xmx{params.mem} -jar {params.jar_file} \
            PE \
            {params.options} \
            -threads {threads} \
            {input.R1} {input.R2} \
            {params.out1P} {params.out1U} {params.out2P} {params.out2U} \
            {params.processing_options} &> {log}

        zcat {params.out1U} {params.out2U} | gzip > {output.out1U} && rm {params.out1U} {params.out2U}
        """

seq_type = "both"

rule map_MT_PE_SE:
    input:
        R1 = "data/reads_filtered/{sample}_{lane}_qc_R1.fastq.gz",
        R2 = "data/reads_filtered/{sample}_{lane}_qc_R2.fastq.gz",
        U = "data/reads_filtered/{sample}_{lane}_qc_U.fastq.gz",
        # R1 = lambda wildcards: read_datasets_inputs("{sample}".format(sample=wildcards.sample), "1"),
        # R2 = lambda wildcards: read_datasets_inputs("{sample}".format(sample=wildcards.sample), "2"),
        gmap_db = gmap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.chromosome"
    output:
        outmt_sam = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        gmap_db = lambda wildcards: wildcards.ref_genome_mt,
        RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample',
        uncompressed_output = lambda wildcards, output: output.outmt_sam.replace("_outmt.sam.gz", "_outmt.sam")
    log:
        log_dir + "/{sample}/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{lane}_{ref_genome_mt}_map_MT_PE_SE.log"
    #conda: "envs/environment.yaml"
    threads:
        config["map"]["gmap_threads"]
    message: "Mapping reads for read dataset {wildcards.sample}_{wildcards.lane} to {wildcards.ref_genome_mt} mt genome"
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
        outmt_sam = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt.sam.gz"
        #outmt_sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_outmt.sam.gz"
    output:
        outmt1 = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt1.fastq.gz",
        outmt2 = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt2.fastq.gz",
        outmt = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt.fastq.gz",
        # outmt1 = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_outmt1.fastq.gz",
        # outmt2 = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_outmt2.fastq.gz",
        # outmt = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_outmt.fastq.gz",
        #log = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/sam2fastq.done"
    #conda: "envs/environment.yaml"
    message:
        "Converting {input.outmt_sam} to FASTQ"
    run:
        sam2fastq(samfile = input.outmt_sam, outmt1 = output.outmt1, outmt2 = output.outmt2, outmt = output.outmt)

rule map_nuclear_MT_SE:
    input:
        outmt = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt.fastq.gz",
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome"
    output:
        outS = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_outS.sam.gz"
        # outS = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_outS.sam.gz"
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gmap_db = lambda wildcards, input: os.path.split(input.gmap_db)[1].replace(".chromosome", ""),
        uncompressed_output = lambda wildcards, output: output.outS.replace("_outS.sam.gz", "_outS.sam")
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
        logS = log_dir + "/{sample}/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_map_nuclear_MT_SE.log"
    message:
        "Mapping onto complete human genome (nuclear + mt)... SE reads"
    run:
        if os.path.isfile(input.outmt):
            #"gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} | gzip -c - > {output.outmt_sam} 2> {log}"
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt} &> {log.logS} && gzip {params.uncompressed_output} &>> {log.logS}")
        else:
            open(output.outS, 'a').close()

rule map_nuclear_MT_PE:
    input:
        outmt1 = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt1.fastq.gz",
        outmt2 = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt2.fastq.gz",
        gmap_db = gmap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.chromosome"
    output:
        outP = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_outP.sam.gz"
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
        logP = log_dir + "/{sample}/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_map_nuclear_MT_PE.log"
    message:
        "Mapping onto complete human genome (nuclear + mt)... PE reads"
    run:
        if os.path.isfile(input.outmt1):
            shell("gsnap -D {params.gmap_db_dir} -d {params.gmap_db} -o {params.uncompressed_output} --gunzip -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt1} {input.outmt2} &> {log.logP} && gzip {params.uncompressed_output} &>> {log.logP}")
        else:
            open(output.outP, 'a').close()

rule filtering_mt_alignments:
    input:
        outmt = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_outmt.sam.gz",
        outS = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_outS.sam.gz",
        outP = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_outP.sam.gz"
    output:
        sam = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz"
        #sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.sam"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
        #outdir = lambda wildcards, output: os.path.split(output.sam)[0]
    #conda: "envs/environment.yaml"
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
        sam = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz",
    output:
        "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT.bam",
    message: "Converting {input.sam} to {output}"
    log: log_dir + "/{sample}/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/sam2bam.log"
    #group: "variant_calling"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        zcat {input.sam} | samtools view -b -o {output} - &> {log}
        """

def check_tmp_dir(dir):
    if os.getenv("TMP"):
        TMP = os.getenv("TMP")
    else:
        TMP = dir
    return TMP

rule sort_bam:
    input:
        bam = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT.bam"
    output:
        sorted_bam = "results/{sample}/map/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam"
    message: "Sorting {input.bam} to {output.sorted_bam}"
    params:
        TMP = check_tmp_dir(config["tmp_dir"])
    log: log_dir + "/{sample}/OUT_{sample}_{lane}_{ref_genome_mt}_{ref_genome_n}/map/sort_bam.log"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
    shell:
        """
        samtools sort -o {output.sorted_bam} -T {params.TMP} {input.bam} &> {log}
        # samtools sort -o {output.sorted_bam} -T ${{TMP}} {input.bam}
        """

rule merge_bam:
    input:
        sorted_bams = lambda wildcards: get_sample_bamfiles(datasets_tab, res_dir="results", sample = wildcards.sample, ref_genome_mt = wildcards.ref_genome_mt, ref_genome_n = wildcards.ref_genome_n)
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

# rule index_merged_bam:
#     input:
#         merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam"
#     output:
#         merged_bam_index = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam.bai"
#     log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_index_merge_bam.log"
#     message: "Indexing {input.merged_bam}"
#     #conda: "envs/samtools_biopython.yaml"
#     shell:
#         """
#         samtools index {input} {output} &> {log}
#         """

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
    shell:
        """
        java -jar modules/picard.jar CreateSequenceDictionary R={input.mt_n_fasta} O={output.genome_dict}
        """

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
    message: "Realigning indels in {input.merged_bam} with GATK 3.8 - LeftAlignIndels"
    shell:
        """
        java -Xmx6G -jar modules/GenomeAnalysisTK.jar \
            -R {input.mt_n_fasta} \
            -T LeftAlignIndels \
            -I {input.merged_bam} \
            -o {output.merged_bam_left_realigned} \
            --filter_reads_with_N_cigar
        """

rule bam2pileup:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam",
        # sorted_bam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam",
        genome_index = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
    output:
        pileup = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
        # pileup = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
    params:
        genome_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    message: "Generating pileup {output.pileup} from {input.merged_bam}"
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_bam2pileup.log"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
    shell:
        """
        samtools mpileup -B -f {params.genome_fasta} -o {output.pileup} {input.merged_bam} &> {log}
        """

rule pileup2mt_table:
    input:
        pileup = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
        # pileup = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.pileup"
    output:
        mt_table = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-mt_table.txt"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file"))
    message: "Generating mt_table {output.mt_table} from {input.pileup}, ref mt: {params.ref_mt_fasta}"
    #conda: "envs/environment.yaml"
    #group: "variant_calling"
    run:
        pileup2mt_table(pileup=input.pileup, fasta=params.ref_mt_fasta, mt_table=output.mt_table)

rule make_single_VCF:
    input:
        merged_bam = "results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam",
        # sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz",
        mt_table = "results/{sample}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-mt_table.txt"
        # sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT.sam.gz",
        # mt_table = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/variant_calling/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-mt_table.txt"
    output:
        single_vcf = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz",
        single_bed = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.bed"
    params:
        ref_mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file = get_mt_fasta(reference_tab, wildcards.ref_genome_mt, "ref_genome_mt_file")),
        TMP = check_tmp_dir(config["tmp_dir"])
    message: "Processing {input.merged_bam} to get VCF {output.single_vcf}"
    #conda: "envs/samtools_biopython.yaml"
    #group: "variant_calling"
    run:
        # function (and related ones) from mtVariantCaller
        # vcf_dict = mtvcf_main_analysis(sam_file = input.sam, mtable_file = input.mt_table, name2 = wildcards.sample)
        tmp_sam = os.path.split(input.merged_bam)[1].replace(".bam", ".sam")
        shell("samtools view {merged_bam} > {tmp_dir}/{tmp_sam}".format(merged_bam = input.merged_bam, tmp_dir = params.TMP, tmp_sam = tmp_sam))
        vcf_dict = mtvcf_main_analysis(sam_file = "{tmp_dir}/{tmp_sam}".format(tmp_dir = params.TMP, tmp_sam = tmp_sam), mtable_file = input.mt_table, name2 = wildcards.sample)
        # ref_genome_mt will be used in the VCF descriptive field
        # seq_name in the VCF data
        seq_name = get_seq_name(params.ref_mt_fasta)
        VCF_RECORDS = VCFoutput(vcf_dict, reference = wildcards.ref_genome_mt, seq_name = seq_name, vcffile = output.single_vcf)
        BEDoutput(VCF_RECORDS, seq_name = seq_name, bedfile = output.single_bed)

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
        single_vcf_list = lambda wildcards: get_genome_single_vcf_files(analysis_tab, ref_genome_mt = wildcards.ref_genome_mt),
        index_vcf = lambda wildcards: get_genome_single_vcf_index_files(analysis_tab, ref_genome_mt = wildcards.ref_genome_mt),
        #single_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz",
        #index_vcf = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz.csi"
    output:
        merged_vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.vcf"
    #message: lambda wildcards: "Merging vcf files for mt reference genome: {ref_genome_mt}".format(ref_genome_mt = wildcards.ref_genome_mt)
    message: "Merging vcf files for mt reference genome: {wildcards.ref_genome_mt}"
    #conda: "envs/bcftools.yaml"
    run:
        shell("bcftools merge {input.single_vcf_list} -O v -o {output.merged_vcf}")
