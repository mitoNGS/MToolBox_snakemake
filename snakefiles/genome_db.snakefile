import os, sys
import snakemake
from types import SimpleNamespace
import snakemake, os, sys
from pathlib import Path
import yaml
#import wget
import tarfile

#configfile: "config.yaml"

for path in sys.path:
    if "snakefiles" in path:
        sys.path.append(path.replace("/snakefiles", ""))
        rootdir = path.replace("/snakefiles", "")

from modules.genome_db import (
    check_ref_organism, get_gmap_build_nuclear_mt_input,
    make_ref_organism_dict, run_gmap_build
)

from modules.config_parsers import (
    get_analysis_species, parse_config_tabs, parse_config_tab, get_genome_files
)

genome_db_data_file = os.path.join(rootdir, "data/genome_dbs.yaml")
#print(genome_db_data_file)
reference_tab_file =  os.path.join(rootdir, "data/reference_genomes.tab")

with open(genome_db_data_file) as file:
    genome_db_data = yaml.full_load(file)

reference_tab = parse_config_tab(tab_file=reference_tab_file, index=["ref_organism"])
genome_fasta_dir = os.path.join(rootdir, "data/genomes")

# not sure it should be like that
if config == False:
    # this is not run as subworkflow
    configfile: "config.yaml"
    # res_dir = config["results"]
    # map_dir = config["map_dir"]
    if config["ref_organism"]:
        ref_organism_config = config["ref_organism"].split(",")
    else:
        sys.exit("Please provide at least one ref_organism.")
else:
    # this is a part of a workflow, likely variant_calling
    analysis_tab = parse_config_tab(tab_file="data/analysis.tab", index=["sample"])
    ref_organism_config, analysis_tab = check_ref_organism(config=config, analysis_tab=analysis_tab, reference_tab=reference_tab)

gmap_db_dir = os.path.join(rootdir, config["map"]["gmap_db_dir"])
log_dir = config["log_dir"]

# Build ref_organism_dict. This will be the source
# for final output files of the pipeline.
ref_organism_dict = make_ref_organism_dict(ref_organism_config=ref_organism_config,
                                            gmap_db_dir=gmap_db_dir, genome_db_data=genome_db_data,
                                            genome_fasta_dir=genome_fasta_dir, reference_tab=reference_tab)

#print(ref_organism_dict)

## toy example of ref_organism_dict
## status MUST be "download" or "new"
# ref_organism_dict = {"ggallus" : SimpleNamespace(ref_genome_mt="NC_001224.1",
#                                                     ref_genome_n="GCF_000146045.2_R64_genomic",
#                                                     ref_genome_mt_file="NC_001224.1.fa",
#                                                     ref_genome_n_file="GCF_000146045.2_R64_genomic.fa",
#                                                     species="ggallus",
#                                                     status="download",
#                                                     zenodo_id="bwe")}

def get_gmap_db_outputs(ref_organism_dict=None):
    gmap_db_mt_outputs = []
    for ref_organism in ref_organism_dict:
        gmap_db_mt_outputs.append(gmap_db_dir + "/{ref_organism}/{ref_organism}.done".format(ref_organism=ref_organism))
    return gmap_db_mt_outputs

def get_genome_indexes_dicts_outputs(ref_organism_dict=None):
    genome_indexes_dicts_outputs = []
    for ref_organism in ref_organism_dict:
        genome_indexes_dicts_outputs.append(rootdir + "/data/genomes/{ref_organism}.dict".format(ref_organism=ref_organism))
        genome_indexes_dicts_outputs.append(rootdir + "/data/genomes/{ref_organism}_mt_n.fasta.fai".format(ref_organism=ref_organism))
    return genome_indexes_dicts_outputs
    

# a target rule to define the desired final output
rule all:
    input:
        get_gmap_db_outputs(ref_organism_dict),
        get_genome_indexes_dicts_outputs(ref_organism_dict)

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint check_gmap_db_status:
    input:
        ancient(gmap_db_dir + "/{ref_organism}")
    output:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
        # fetch_mt = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetch",
        # fetch_n = gmap_db_dir + "/{ref_organism}/{ref_organism}_n.fetch"
    params:
        status = lambda wildcards: ref_organism_dict[wildcards.ref_organism].status,
        fetch_mt = lambda wildcards: ref_organism_dict[wildcards.ref_organism].fetch_mt_genome,
        fetch_n = lambda wildcards: ref_organism_dict[wildcards.ref_organism].fetch_n_genome
    run:
        # simulate some output value
        shell("echo 'status: {params.status}' > {output.ref_organism_flag}")
        shell("echo 'fetch_mt_genome: {params.fetch_mt}' >> {output.ref_organism_flag}")
        shell("echo 'fetch_n_genome: {params.fetch_n}' >> {output.ref_organism_flag}")        # shell("echo {params.fetch_n} > {output.fetch_n}")
        #"echo {params.status} > {output.mt_n}"
        # create genome "status" file

# ## wanna dl genomes?
# def evaluate_genome_fasta_mt(wildcards):
#     with open(checkpoints.check_gmap_db_status.get(ref_organism=wildcards.ref_organism).output.fetch_mt) as f:
#         if f.read().strip() == "yes":
#             return gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetched", \
#                 gmap_db_dir + "/{ref_organism}/{ref_organism}_mt_n.fetched"
#         else:
#             return gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.new", \
#                 gmap_db_dir + "/{ref_organism}/{ref_organism}_mt_n.new"

rule download_mt_genome_fasta:
    input:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
        #ref_organism_fetch = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetch",
    output:
        ref_organism_fetch = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetched",
#        ref_organism_fetch = "data/genomes/{ref_genome_mt_file}",
    params:
        mt_genome_fasta    = lambda wildcards: "{rootdir}/data/genomes/{ref_genome_mt_file}".format(rootdir=rootdir,
                                                                    ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file),
        ref_genome_mt      = lambda wildcards: "{ref_genome_mt}".format(ref_genome_mt=ref_organism_dict[wildcards.ref_organism].ref_genome_mt),
        #ref_genome_mt_file = lambda wildcards: "{ref_genome_mt_file}".format(ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file)
    run:
        gmap_db_status = yaml.full_load(open(input.ref_organism_flag, 'r'))
        os.makedirs("{genome_fasta_dir}".format(genome_fasta_dir=os.path.split(params.mt_genome_fasta)[0]), exist_ok=True)
        #print(input.ref_organism_flag, gmap_db_status)
        if gmap_db_status["fetch_mt_genome"] == True:
            shell("esearch -db nucleotide -query '{params.ref_genome_mt}' | efetch -format fasta > {params.mt_genome_fasta} && touch {output.ref_organism_fetch}")
        else:
            shell("touch {output.ref_organism_fetch}")
        #shell("echo done > {output.ref_organism_fetch}")

rule download_mt_n_genome_fasta:
    input:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
        #ref_organism_fetch = gmap_db_dir + "/{ref_organism}/{ref_organism}_n.fetch",
    output:
        ref_organism_fetch = gmap_db_dir + "/{ref_organism}/{ref_organism}_n.fetched",
        #ref_organism_fetch = "data/genomes/{ref_genome_n_file}",
    params:
        n_genome_fasta = lambda wildcards: "{rootdir}/data/genomes/{ref_genome_n_file}".format(rootdir=rootdir,
                                                                    ref_genome_n_file=ref_organism_dict[wildcards.ref_organism].ref_genome_n_file),
        ref_genome_n = lambda wildcards: "{ref_genome_n}".format(ref_genome_n=ref_organism_dict[wildcards.ref_organism].ref_genome_n),
        #ref_genome_n_file = lambda wildcards: "{ref_genome_mt}".format(ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_n_file)
    run:
        gmap_db_status = yaml.full_load(open(input.ref_organism_flag, 'r'))
        os.makedirs("{genome_fasta_dir}".format(genome_fasta_dir=os.path.split(params.n_genome_fasta)[0]), exist_ok=True)
        if gmap_db_status["fetch_n_genome"] == True:
            if ref_organism_dict[wildcards.ref_organism].ref_genome_n.startswith("GCA"):
                db = "insdc"
            elif ref_organism_dict[wildcards.ref_organism].ref_genome_n.startswith("GCF"):
                db = "refseq"
            else:
                sys.exit("Please provide a valid assembly accession (GCF- or GCA-).")
            # this
            shell("esearch -db assembly -query '{ref_genome_n}' | elink -target nucleotide -name assembly_nuccore_{db} | efetch -format fasta > {n_genome_fasta} && touch {ref_organism_fetch}".format(ref_genome_n=params.ref_genome_n, db=db, n_genome_fasta=params.n_genome_fasta, ref_organism_fetch=output.ref_organism_fetch))
            # or this
            #shell("esearch -db nucleotide -query '{params.ref_genome_n}' | efetch -format fasta > {params.mt_genome_fasta} && touch {output.ref_organism_fetch}")
        else:
            shell("touch {output.ref_organism_fetch}")

# workflow, option1
rule make_mt_gmap_db:
    input:
#        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
        ref_organism_fetch_mt = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetched",
        #mt_genome_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file)
    output:
        dl_mt = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.new",
    params:
        gmap_db_dir = lambda wildcards, output: os.path.split(output.dl_mt)[0],
        gmap_db = lambda wildcards: "{ref_genome_mt}".format(ref_genome_mt=ref_organism_dict[wildcards.ref_organism].ref_genome_mt),
        mt_genome_fasta = lambda wildcards: "{rootdir}/data/genomes/{ref_genome_mt_file}".format(rootdir=rootdir,
                                                                    ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file)
        #gmap_db = lambda wildcards, output: os.path.split(output.gmap_db)[1].replace(".chromosome", "")
    message: "Generating gmap db for mt genome: {params.mt_genome_fasta}.\nWildcards: {wildcards}"
    log: "logs/gmap_build/{ref_organism}_mt.log"
    #conda: "envs/environment.yaml"
    run:
        run_gmap_build(mt_genome_file=params.mt_genome_fasta, gmap_db_dir=params.gmap_db_dir,
                        gmap_db=params.gmap_db, log=log, mt_is_circular=True)
        shell("touch {output.dl_mt}")

rule get_gmap_build_nuclear_mt_input:
    input:
        ref_organism_fetch_mt = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.fetched",
        ref_organism_fetch_n = gmap_db_dir + "/{ref_organism}/{ref_organism}_n.fetched",
        # ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
        # mt_genome_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file),
        # n_genome_fasta = lambda wildcards: "data/genomes/{ref_genome_n_file}".format(ref_genome_n_file=ref_organism_dict[wildcards.ref_organism].ref_genome_n_file),
    output:
        mt_n_fasta = rootdir + "/data/genomes/{ref_organism}_mt_n.fasta"
    params:
        n_genome_fasta = lambda wildcards: "{rootdir}/data/genomes/{ref_genome_n_file}".format(rootdir=rootdir,
                                                                    ref_genome_n_file=ref_organism_dict[wildcards.ref_organism].ref_genome_n_file),
        mt_genome_fasta = lambda wildcards: "{rootdir}/data/genomes/{ref_genome_mt_file}".format(rootdir=rootdir,
                                                                    ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file),
    run:
        get_gmap_build_nuclear_mt_input(n_genome_file=params.n_genome_fasta, mt_genome_file=params.mt_genome_fasta, n_mt_file=output.mt_n_fasta)


rule make_mt_n_gmap_db:
    input:
        mt_n_fasta = rootdir + "/data/genomes/{ref_organism}_mt_n.fasta",
        #mt_fasta = lambda wildcards: "data/genomes/{ref_genome_mt_file}".format(ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file),
        # mt_n_fasta = rules.get_gmap_build_nuclear_mt_input.output.mt_n_fasta,
    output:
        dl_mt_n = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt_n.new",
    params:
        gmap_db_dir = lambda wildcards, output: os.path.split(output.dl_mt_n)[0],
        gmap_db = lambda wildcards: "{ref_genome_mt}_{ref_genome_n}".format(ref_genome_mt=ref_organism_dict[wildcards.ref_organism].ref_genome_mt,
                                                                            ref_genome_n=ref_organism_dict[wildcards.ref_organism].ref_genome_n),
        mt_fasta = lambda wildcards: "{rootdir}/data/genomes/{ref_genome_mt_file}".format(rootdir=rootdir,
                                                                    ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file)
    message: "Generating joint gmap db for mt+n genome: {wildcards.ref_organism}"
    log: "logs/gmap_build/{ref_organism}_mt_n.log"
    run:
        run_gmap_build(mt_n_genome_file=input.mt_n_fasta, mt_genome_file=params.mt_fasta,
                            gmap_db_dir=params.gmap_db_dir, gmap_db=params.gmap_db, log=log, mt_is_circular=True)
        shell("touch {output}")

# workflow, option2
rule download_genome_db:
    input:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
    output:
        dl_mt = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.download",
        dl_mt_n = gmap_db_dir + "/{ref_organism}/{ref_organism}_mt_n.download",
    params:
        zenodo_id = lambda wildcards: genome_db_data[wildcards.ref_organism]["zenodo_id"],
    message:
        "Downloading genome"
    run:
        # download archive
        shell("zenodo_get -w /tmp/{params.zenodo_id}.list {params.zenodo_id}")
        dl_filename = os.path.split(open("/tmp/{}.list".format(params.zenodo_id), 'r').read().strip())[1]
        shell("zenodo_get -o /tmp {params.zenodo_id}")
        # extract it
        tar = tarfile.open("/tmp/{}".format(dl_filename))
        tar.extractall(path="{gmap_db_dir}/{ref_organism}".format(gmap_db_dir=gmap_db_dir, ref_organism=wildcards.ref_organism))
        tar.close()
        #shell("touch {output.dl_mt}")
        shell("touch {output.dl_mt_n}")
        shell("touch {output.dl_mt}")

# input function for the rule aggregate
def aggregate_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with open(checkpoints.check_gmap_db_status.get(ref_organism=wildcards.ref_organism).output.ref_organism_flag) as f:
        gmap_db_status = yaml.full_load(f)
        if gmap_db_status["status"] == "download":
#        if f.read().strip() == "download":
            return gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.download", \
                gmap_db_dir + "/{ref_organism}/{ref_organism}_mt_n.download"
        else:
            return gmap_db_dir + "/{ref_organism}/{ref_organism}_mt.new", \
                gmap_db_dir + "/{ref_organism}/{ref_organism}_mt_n.new"

rule index_genome:
    input:
        mt_n_fasta = rules.get_gmap_build_nuclear_mt_input.output.mt_n_fasta
        #mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        # genome_index = rootdir + "/data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta.fai"
        genome_index = rootdir + "/data/genomes/{ref_organism}_mt_n.fasta.fai"
    message: "Indexing {input.mt_n_fasta} with samtools faidx"
    # log: log_dir + "/{ref_genome_mt}_{ref_genome_n}.samtools_index.log"
    log: log_dir + "/{ref_organism}.samtools_index.log"
    #conda: "envs/samtools_biopython.yaml"
    shell:
        """
        samtools faidx {input.mt_n_fasta} &> {log}
        """

rule dict_genome:
    input:
        mt_n_fasta = rules.get_gmap_build_nuclear_mt_input.output.mt_n_fasta
        #mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        # genome_dict = rootdir + "/data/genomes/{ref_genome_mt}_{ref_genome_n}.dict"
        genome_dict = rootdir + "/data/genomes/{ref_organism}.dict"
    message: "Creating .dict of {input.mt_n_fasta} with picard CreateSequenceDictionary"
    log: log_dir + "/{ref_organism}.picard_dict.log"
    #conda: "envs/samtools_biopython.yaml"
    run:
        shell("picard CreateSequenceDictionary R={input.mt_n_fasta} O={output.genome_dict} &> {log}")

rule check_gmap_db_status_2: # aggregate
    input:
        aggregate_input
    output:
        mt = gmap_db_dir + "/{ref_organism}/{ref_organism}.done",
#        mt_n = gmap_db_dir + "{ref_organism}/{ref_organism}.done"
    shell:
        "touch {output.mt}"

