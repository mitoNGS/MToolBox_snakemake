#!/usr/bin/env python
import gzip
import os
import time

import pandas as pd
from sqlalchemy import create_engine
from types import SimpleNamespace

from modules.general import get_SAM_header, memory_usage_resource


# TODO: ext and n_occurrences are not used anywhere
def read_sam_file_only_readID_chunks_intoSQL(samfile,
                                             n_occurrences=1,
                                             chunksize=100000,
                                             table_name="outS",
                                             ext=".sam.gz"):
    """ Read a SAM file, then keep a list of IDs of reads occurring
        <n_occurrences> times.
        In the specific case of outS.sam and outP.sam files, these are the IDs
        of the reads we want to keep in the OUT.sam file.
        Load the entries in a SQL db.
    """
    # TODO: n and n_occurrences are not used anywhere
    n = n_occurrences

    # TODO: in-memory db is really good?!
    # Create in-memory SQLite db
    engine = create_engine('sqlite://', echo=False)
    # samfile = path/to/out.sam --> table_name = out
    # TODO: what does this mean?
    # table_name = table_name

    # Read the SAM file in chunks
    header_lines, comment_count = get_SAM_header(samfile)
    # function that reads a samfile and skips rows with unaligned reads
    t = pd.read_table(samfile,
                      sep='\t',
                      skiprows=comment_count,
                      chunksize=chunksize,
                      usecols=[0, 2],
                      names=['readID', 'RNAME'],
                      compression="infer",
                      index_col=False)

    for chunk in t:
        elapsed = time.time()
        chunk = chunk.query('RNAME != "*"')
        chunk = chunk.drop(columns=['RNAME'])
        chunk.to_sql(table_name, con=engine, if_exists="append")
        print("{} seconds, memory: {} MB".format(time.time() - elapsed,
                                                 memory_usage_resource()))

    return engine, table_name


# TODO: ref_mt_fasta is not used anywhere
def filter_alignments(outmt=None, outS=None, outP=None, OUT=None,
                      ref_mt_fasta=None):
    print("Processing {}".format(outS))
    outS_sql, table_name_S = read_sam_file_only_readID_chunks_intoSQL(outS,
                                                                      table_name="outS")
    print("Processing {}".format(outP))
    outP_sql, table_name_P = read_sam_file_only_readID_chunks_intoSQL(outP,
                                                                      table_name="outP")

    good_reads_S = pd.read_sql_query(
        "SELECT readID FROM {} GROUP BY readID HAVING COUNT(*) == 1".format(table_name_S),
        outS_sql, chunksize=100000)
    print("SQL query on outS, memory: {} MB".format(memory_usage_resource()))
    good_reads_P = pd.read_sql_query(
        "SELECT readID FROM {} GROUP BY readID HAVING COUNT(*) == 2".format(table_name_P),
        outP_sql, chunksize=100000)
    print("SQL query on outP, memory: {} MB".format(memory_usage_resource()))
    # TODO: the following is not good for performance
    good_reads = pd.DataFrame()
    for c in good_reads_P:
        good_reads = good_reads.append(c)
    print("good_reads_P append, memory: {} MB".format(memory_usage_resource()))
    for c in good_reads_S:
        good_reads = good_reads.append(c)
    print("good_reads_S append, memory: {} MB".format(memory_usage_resource()))
    print("Total reads to extract alignments of: {}".format(len(good_reads)))

    samfile = outmt
    tc = pd.read_table(samfile,
                       sep='\t',
                       skiprows=get_SAM_header(samfile)[1],
                       chunksize=100000,
                       header=None,
                       engine="python",
                       names=["readID", "FLAG", "RNAME"] + list("QWERTYUIOPASDFGHJK"),
                       compression="infer",
                       index_col=False)

    # open OUT.sam file and write SAM header from outS.sam (outP would be the same).
    OUT_uncompressed = OUT.replace(".gz", "")
    with open(OUT_uncompressed, "w") as f:
        sss = gzip.open(outS, 'rb')
        l = sss.readline().decode("utf-8")
        while l[0] == "@":
            if not l.startswith("@PG"):
                f.write(l)
            l = sss.readline().decode("utf-8")
        f.write("\t".join(["@RG", "ID:sample", "PL:illumina", "SM:sample"]) + "\n")

    n_extracted_alignments = 0
    for chunk in tc:
        chunk = chunk.query('RNAME != "*"')
        OUT_chunk = pd.merge(chunk, good_reads, how="inner", on="readID")
        print("Chunk, memory: {} MB".format(memory_usage_resource()))
        n_extracted_alignments += len(OUT_chunk)
        # Append alignments to OUT.sam
        OUT_chunk.to_csv(OUT_uncompressed, mode="a", header=False, sep="\t",
                         index=False)

    print("Compressing OUT.sam file")
    os.system("gzip {}".format(OUT_uncompressed))
    print("OUT.sam compressed, memory: {} MB".format(memory_usage_resource()))
    print("Total alignments extracted: {}".format(n_extracted_alignments))

def cat_alignment(samfile=None, outfile=None, ref_mt_fasta_header=None):
    #samhandle = gzip.open(samfile, 'rt')
    samhandle = open(samfile, 'r')
    outhandle = gzip.open(outfile, 'at')
    total_alignments = 0
    filtered_alignments = 0
    for l in samhandle:
        total_alignments += 1
        if l.startswith("@") == False and l.split()[2] == ref_mt_fasta_header:
            filtered_alignments += 1
            outhandle.write(l)
    samhandle.close()
    outhandle.close()
    return total_alignments, filtered_alignments

def cat_alignments(*samfiles, outfile=None, ref_mt_fasta_header=None):
    filtering_report = {}
    header_lines, comment_count = get_SAM_header(samfiles[0][0])
    outhandle = gzip.open(outfile, 'at')
    for l in header_lines:
        outhandle.write(l)
    outhandle.close()
    for samfile in samfiles[0]:
        print(samfile)
        total_alignments, filtered_alignments = cat_alignment(samfile=samfile, outfile=outfile, ref_mt_fasta_header=ref_mt_fasta_header)
        filtering_report[samfile] = SimpleNamespace(samfile=samfile, total_alignments=total_alignments, filtered_alignments=filtered_alignments, ref=ref_mt_fasta_header)
    return filtering_report