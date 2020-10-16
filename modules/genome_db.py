#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from modules.general import is_compr_file
from Bio import SeqIO, bgzf
from snakemake import shell
from types import SimpleNamespace
import gzip
import os

shell.prefix("set -euo pipefail;")

def open_genome_file(genome_file=None):
    if is_compr_file(genome_file):
        return SeqIO.parse(gzip.open(genome_file, 'rt'), 'fasta')
    else:
        return SeqIO.parse(genome_file, 'fasta')

def check_mt_in_n(s_n_handle=None, mt_genome_id=None):
    """
    s_n_handle is a SeqRecord instance
    mt_genome_id is a string
    """
    mt_in_n = False
    if "mitochondri" in s_n_handle.description:
        mt_in_n = True
    elif s_n_handle.id == mt_genome_id:
        mt_in_n = True
    elif "chrM" in s_n_handle.id or "chrMT" in s_n_handle.id:
        mt_in_n = True
    return mt_in_n

def get_gmap_build_nuclear_mt_input(n_genome_file=None, mt_genome_file=None, n_mt_file=None):
    """Takes a Bio.SeqIO parsed nuclear genome and mt genome,
    checks if the mt genome is already in the nuclear genome handle,
    generates fasta file with genome for gmap_build.
    
    PLEASE NOTE: atm it's safe to have one-contig mt genomes.
    
    Args:
        n_handle:       fasta file
        mt_handle:      fasta file
    
    Return:
        ###
    """
    n_handle = open_genome_file(genome_file=n_genome_file)
    mt_handle = open_genome_file(genome_file=mt_genome_file)
    mt_n_fasta = open(n_mt_file, 'w')
    #mt_n_fasta = bgzf.BgzfWriter(n_mt_file, 'w')
    for s in mt_handle:
        mt_genome_id = s.id
        mt_genome_seq = str(s.seq)
    for s in n_handle:
        # if s.id != mt_genome_id:
        if check_mt_in_n(s_n_handle=s, mt_genome_id=mt_genome_id) == False:
            mt_n_fasta.write(">{}\n{}\n".format(s.id, str(s.seq)))
    mt_n_fasta.write(">{}\n{}\n".format(mt_genome_id, mt_genome_seq))
    mt_n_fasta.close()
    #return True

def get_mt_header(mt_genome_file=None):
    """Gets seq id of a single-contig fasta file"""
    mt_handle = SeqIO.read(mt_genome_file, 'fasta')
    return mt_handle.id    

def run_gmap_build(mt_n_genome_file=None, mt_genome_file=None,
                    gmap_db_dir=None, gmap_db=None, log=None, mt_is_circular=True):
    """
    gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -g -s none {output.mt_n_fasta} 2> /dev/null | gmap_build -D {params.gmap_db_dir} -d {params.gmap_db} -s none {output.mt_n_fasta} &> {log}
    """
    #print("Input files provided: n_genome_file={}, mt_genome_file={}".format(n_genome_file, mt_genome_file))
    # nuclear + mt db
    c_flag = ""
    g_flag = ""
    if mt_is_circular:
        mt_id = get_mt_header(mt_genome_file=mt_genome_file)
        c_flag = "-o {}".format(mt_id)
    if mt_n_genome_file:
        #get_gmap_build_nuclear_mt_input(n_genome_file=n_genome_file, mt_genome_file=mt_genome_file, n_mt_file=n_mt_file)
        shell("gmap_build -D {gmap_db_dir} -d {gmap_db} {c_flag} -s none {input_fasta} &> {log}".format(gmap_db_dir=gmap_db_dir,
                                                                                            gmap_db=gmap_db, c_flag=c_flag, input_fasta=mt_n_genome_file,
                                                                                            log=log))
    # mt db
    else:
        if is_compr_file(mt_genome_file):
            g_flag = "-g"
        # else:
        #     g_flag = ""
        shell("gmap_build -D {gmap_db_dir} -d {gmap_db} {g_flag} {c_flag} -s none {input_fasta} &> {log}".format(gmap_db_dir=gmap_db_dir,
                                                                                        gmap_db=gmap_db, input_fasta=mt_genome_file,
                                                                                        log=log, g_flag=g_flag, c_flag=c_flag))

# if species is not defined by config.yaml, should be parsed for each analysis
def check_ref_organism(config=None, analysis_tab=None, reference_tab=None):
    ref_organism_config = config["ref_organism"]
    if ref_organism_config is not None:
        # ref_organism is defined through config (either config.yaml or command line)
        if len(ref_organism_config.split(",")) > 1:
            sys.exit("Please provide only one reference organism (--config ref_organism)")
        else:
            ref_organism = ref_organism_config
            if ref_organism in list(reference_tab["ref_organism"]) or ref_organism in genome_db_data:
                analysis_tab = analysis_tab.assign(ref_organism=ref_organism)
            else:
                sys.exit("Provided species {} in not present in reference_genomes.tab.".format(species))
    else:
        # ref_organism is defined through analysis_tab
        ref_organism_config = list(set(analysis_tab['ref_organism']))
        if ref_organism_config is None:
            sys.exit("Please provide at least one reference organism (--config ref_organism)")
    return ref_organism_config, analysis_tab

def make_ref_organism_dict(ref_organism_config=None, gmap_db_dir=None,
                            genome_db_data=None, genome_fasta_dir=None,
                            reference_tab=None):
    ref_organism_dict = {}
    for r in ref_organism_config:
        if os.path.isdir("{gmap_db_dir}/{ref_organism}".format(gmap_db_dir=gmap_db_dir,
                                                            ref_organism=r)):
            pass
        else:
            os.makedirs("{gmap_db_dir}/{ref_organism}".format(gmap_db_dir=gmap_db_dir,
                                                            ref_organism=r), exist_ok=True)
        ref_organism_dict[r] = SimpleNamespace()
        # check if it's in genome_db_data
        if r in genome_db_data:
            os.makedirs("{gmap_db_dir}/{ref_organism}".format(gmap_db_dir=gmap_db_dir,
                                                                        ref_organism=r), exist_ok=True)
            for f in genome_db_data[r]:
                setattr(ref_organism_dict[r], f, genome_db_data[r][f])
                setattr(ref_organism_dict[r], "status", "download")
                setattr(ref_organism_dict[r], "fetch_mt_genome", "no")
                setattr(ref_organism_dict[r], "fetch_n_genome", "no")
                # otherwise in reference_tab
        elif r in reference_tab.index:
            # reference_tab.loc["ggallus_2"]["ref_genome_mt"]
            for attribute in ["ref_genome_mt", "ref_genome_n"]:
                # these attributes MUST be set
                try:
                    a = reference_tab.loc[r][attribute]
                except KeyError:
                    sys.exit("{r} doesn't have a valid {attribute}".format(r=r, attribute=attribute))
                setattr(ref_organism_dict[r], attribute, a)
            for attribute in ["ref_genome_mt_file", "ref_genome_n_file"]:
                # these attributes might not be set and their values
                # will be set based on their related attributes, ie
                # ref_genome_mt_file from ref_genome_mt
                # ref_genome_n_file  from ref_genome_n
                try:
                    a = reference_tab.loc[r][attribute]
                except KeyError:
                    # eg if ref_genome_mt_file is not set in the reference_tab,
                    # its value will be set as <ref_genome_mt>.fasta
                    a = "{ref}.fasta".format(ref=reference_tab.loc[r][attribute.replace("_file", "")])
                setattr(ref_organism_dict[r], attribute, a)
            setattr(ref_organism_dict[r], "status", "new")
            # for attribute in ["ref_genome_mt", "ref_genome_n", "ref_genome_mt_file", "ref_genome_n_file"]:
            #     try:
            #         a = reference_tab.loc[r][attribute]
            #     except KeyError:
            #         sys.exit("{r} doesn't have a valid {attribute}".format(r=r, attribute=attribute))
            #     setattr(ref_organism_dict[r], attribute, a)
            # setattr(ref_organism_dict[r], "status", "new")
            if os.path.isfile(os.path.join(genome_fasta_dir, ref_organism_dict[r].ref_genome_mt_file)):
                setattr(ref_organism_dict[r], "fetch_mt_genome", "no")
                #open(gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetched".format(ref_organism=r), 'a').close()
            else:
                setattr(ref_organism_dict[r], "fetch_mt_genome", "yes")
            if os.path.isfile(os.path.join(genome_fasta_dir, ref_organism_dict[r].ref_genome_n_file)):
                setattr(ref_organism_dict[r], "fetch_n_genome", "no")
                #open(gmap_db_dir + "/{ref_organism}/{ref_organism}_n.fetched".format(ref_organism=r), 'a').close()
            else:
                setattr(ref_organism_dict[r], "fetch_n_genome", "yes")
        # otherwise drop it
        else:
            "{r} is not a reference organism neither in the default genome_db nor in the reference_genomes.tab file. It will be discarded.".format(r=r)
            pass
    return ref_organism_dict
