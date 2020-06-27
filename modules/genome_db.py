#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from modules.general import is_compr_file
from Bio import SeqIO, bgzf
from snakemake import shell
import gzip

shell.prefix("set -euo pipefail;")

def open_genome_file(genome_file=None):
    if is_compr_file(genome_file):
        return SeqIO.parse(gzip.open(genome_file, 'rt'), 'fasta')
    else:
        return SeqIO.parse(genome_file, 'fasta')

def get_gmap_build_nuclear_mt_input(n_genome_file=None, mt_genome_file=None, n_mt_file=None):
    """Takes a Bio.SeqIO parsed nuclear genome and mt genome,
    checks if the mt genome is already in the nuclear genome handle,
    generates compressed fasta file with genome for gmap_build.
    
    PLEASE NOTE: atm it's safe to have one-contig mt genomes.
    
    Args:
        n_handle:       fasta file
        mt_handle:      fasta file
    
    Return:
        ###
    """
    n_handle = open_genome_file(genome_file=n_genome_file)
    mt_handle = open_genome_file(genome_file=mt_genome_file)
    mt_n_fasta = bgzf.BgzfWriter(n_mt_file, 'w')
    for s in mt_handle:
        mt_genome_id = s.id
        mt_genome_seq = str(s.seq)
    for s in n_handle:
        if s.id != mt_genome_id:
            mt_n_fasta.write(">{}\n{}\n".format(s.id, str(s.seq)))
    mt_n_fasta.write(">{}\n{}\n".format(mt_genome_id, mt_genome_seq))
    mt_n_fasta.close()
    #return True

def run_gmap_build(n_genome_file=None, mt_genome_file=None, n_mt_file=None,
                    gmap_db_dir=None, gmap_db=None, log=None):
    """
    gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -g -s none {output.mt_n_fasta} 2> /dev/null | gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -s none {output.mt_n_fasta} &> {log}
    """
    #print("Input files provided: n_genome_file={}, mt_genome_file={}".format(n_genome_file, mt_genome_file))
    # nuclear + mt db
    if n_genome_file:
        get_gmap_build_nuclear_mt_input(n_genome_file=n_genome_file, mt_genome_file=mt_genome_file, n_mt_file=n_mt_file)
        shell("gmap_build -D {gmap_db_dir} -d {gmap_db} -g -s none {input_fasta} &> {log}".format(gmap_db_dir=gmap_db_dir,
                                                                                            gmap_db=gmap_db, input_fasta=n_mt_file,
                                                                                            log=log))
    # mt db
    else:
        if is_compr_file(mt_genome_file):
            g_flag = "-g"
        else:
            g_flag = ""
        shell("gmap_build -D {gmap_db_dir} -d {gmap_db} {g_flag} -s none {input_fasta} &> {log}".format(gmap_db_dir=gmap_db_dir,
                                                                                        gmap_db=gmap_db, input_fasta=mt_genome_file,
                                                                                        log=log, g_flag=g_flag))
