#!/usr/bin/env python

import pandas as pd
from sqlalchemy import create_engine
from modules.general import get_SAM_header, memory_usage_resource
import time, os, gzip

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
