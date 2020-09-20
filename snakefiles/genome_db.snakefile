import os, sys
import snakemake
from types import SimpleNamespace
import snakemake, os, sys
from pathlib import Path
import yaml
#import wget
import tarfile

for path in sys.path:
    if "snakefiles" in path:
        sys.path.append(path.replace("/snakefiles", ""))

from modules.config_parsers import (
    get_analysis_species, parse_config_tabs, parse_config_tab, get_genome_files
)

for path in sys.path:
    if "snakefiles" in path:
        rootdir = path.replace("/snakefiles", "/")
        #print(rootdir)

genome_db_data_file = os.path.join(rootdir, "data/genome_dbs.yaml")
#print(genome_db_data_file)
reference_tab_file =  os.path.join(rootdir, "data/reference_genomes.tab")

with open(genome_db_data_file) as file:
    genome_db_data = yaml.full_load(file)

reference_tab = parse_config_tab(tab_file=reference_tab_file, index=["ref_organism"])

print(genome_db_data)
print(genome_db_data["ggallus"]["zenodo_id"])

configfile: "config.yaml"
# res_dir = config["results"]
# map_dir = config["map_dir"]
log_dir = config["log_dir"]
gmap_db_dir = config["map"]["gmap_db_dir"]
ref_organisms_config = config["ref_organism"].split(",")

# Build ref_organism_dict. This will be the source
# for final output files of the pipeline.
ref_organism_dict = {}
for r in ref_organisms_config:
    # check if it's in genome_db_data
    if r in genome_db_data:
        ref_organism_dict[r] = SimpleNamespace()
        for f in genome_db_data[r]:
            setattr(ref_organism_dict[r], f, genome_db_data[r][f])
            setattr(ref_organism_dict[r], "status", "download")
    # otherwise in reference_tab
    elif r in reference_tab.index:
        # reference_tab.loc["ggallus_2"]["ref_genome_mt"]
        for attribute in ["ref_genome_mt", "ref_genome_n", "ref_genome_mt_file", "ref_genome_n_file"]:
            try:
                a = reference_tab.loc[r][attribute]
            except KeyError:
                sys.exit("{r} doesn't have a valid {attribute}".format(r=r, attribute=attribute))
            setattr(ref_organism_dict[r], attribute, a)
        setattr(ref_organism_dict[r], "status", "new")
    # otherwise drop it
    else:
        "{r} is not a reference organism neither in the default genome_db nor in the reference_genomes.tab file. It will be discarded.".format(r=r)
        pass
        
## toy example of ref_organism_dict
## status MUST be "download" or "new"
# ref_organism_dict = {"ggallus" : SimpleNamespace(ref_genome_mt="NC_001224.1",
#                                                     ref_genome_n="GCF_000146045.2_R64_genomic",
#                                                     ref_genome_mt_file="NC_001224.1.fa",
#                                                     ref_genome_n_file="GCF_000146045.2_R64_genomic.fa",
#                                                     species="ggallus",
#                                                     status="download",
#                                                     zenodo_id="bwe")}

for ref_organism in ref_organism_dict:
    os.makedirs("{gmap_db_dir}/{ref_organism}".format(gmap_db_dir=gmap_db_dir,
                                                        ref_organism=ref_organism), exist_ok=True)
    # os.makedirs("{gmap_db_dir}/{ref_organism}/{ref_organism}.status".format(gmap_db_dir=gmap_db_dir,
    #                                                                         ref_organism=ref_organism), exist_ok=True)

def get_gmap_db_outputs(ref_organism_dict=None):
    for ref_organism in ref_organism_dict:
        gmap_db_mt_outputs = []
        gmap_db_mt_n_outputs = []
        gmap_db_mt_outputs.append(gmap_db_dir + "/{ref_organism}/{ref_organism}.done".format(ref_organism=ref_organism))
    return gmap_db_mt_outputs + gmap_db_mt_n_outputs

# a target rule to define the desired final output
rule all:
    input:
        get_gmap_db_outputs(ref_organism_dict)

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint check_gmap_db_status:
    input:
        ref_organism = gmap_db_dir + "/{ref_organism}"
    output:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
    params:
        status = lambda wildcards: ref_organism_dict[wildcards.ref_organism].status
    run:
        # simulate some output value
        shell("echo {params.status} > {output.ref_organism_flag}")
        #"echo {params.status} > {output.mt_n}"
        # create genome "status" file

# workflow, option1
rule make_mt_gmap_db:
    input:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",#.format(ref_organism=wildcards.ref_organism),
#                                                                                            ref_genome_mt=ref_organism_dict[wildcards.ref_organism].ref_genome_mt),
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}", ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file)
    output:
        dl_mt = gmap_db_dir + "/{ref_organism}/{ref_organism}.new",
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
    message: "Generating gmap db for mt genome: {input.mt_genome_fasta}.\nWildcards: {wildcards}"
    log: "logs/gmap_build/{ref_organism}.log"
    #conda: "envs/environment.yaml"
    run:
        run_gmap_build(mt_genome_file=input.mt_genome_fasta, gmap_db_dir=params.gmap_db_dir,
                        gmap_db=params.gmap_db, log=log, mt_is_circular=True)
        shell("touch {output.dl_mt}")


rule get_gmap_build_nuclear_mt_input:
    input:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}",
                                                    ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file),
        n_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_n_file}",
                                                  ref_genome_n_file=ref_organism_dict[wildcards.ref_organism].ref_genome_n_file),
    output:
        mt_n_fasta = "data/genomes/{ref_genome_mt}_{ref_genome_n}.fasta"
    run:
        get_gmap_build_nuclear_mt_input(n_genome_file=input.n_genome_fasta, mt_genome_file=input.mt_genome_fasta, n_mt_file=output.mt_n_fasta)


rule make_mt_n_gmap_db:
    input:
        mt_genome_fasta = lambda wildcards: expand("data/genomes/{ref_genome_mt_file}", ref_genome_mt_file=ref_organism_dict[wildcards.ref_organism].ref_genome_mt_file),
        mt_n_fasta = rules.get_gmap_build_nuclear_mt_input.output.mt_n_fasta,
    output:
        dl_mt_n = gmap_db_dir + "/{ref_organism}/{ref_organism}.new",
    params:
        gmap_db_dir = config["map"]["gmap_db_dir"],
    message: "Generating gmap db for mt + n genome: {ref_organism}"
    log: "logs/gmap_build/{ref_organism}.log"
    run:
        run_gmap_build(mt_n_genome_file=input.mt_n_fasta, mt_genome_file=input.mt_fasta,
                            gmap_db_dir=params.gmap_db_dir, gmap_db=params.gmap_db, log=log, mt_is_circular=True)
        shell("touch {output}")

# workflow, option2
rule download_genome_db:
    input:
        ref_organism_flag = gmap_db_dir + "/{ref_organism}/{ref_organism}.status",
    output:
        dl_mt_n = gmap_db_dir + "/{ref_organism}/{ref_organism}.download",
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
        tar.extractall(path=gmap_db_dir + wildcards.ref_organism)
        tar.close()
        #shell("touch {output.dl_mt}")
        shell("touch {output.dl_mt_n}")

# input function for the rule aggregate
def aggregate_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with open(checkpoints.check_gmap_db_status.get(ref_organism=wildcards.ref_organism).output.ref_organism_flag) as f:
        if f.read().strip() == "download":
            return gmap_db_dir + "/{ref_organism}/{ref_organism}.download", \
                gmap_db_dir + "/{ref_organism}/{ref_organism}.download"
        else:
            return gmap_db_dir + "/{ref_organism}/{ref_organism}.new", \
                gmap_db_dir + "/{ref_organism}/{ref_organism}.new"

rule check_gmap_db_status_2: # aggregate
    input:
        aggregate_input
    output:
        mt = gmap_db_dir + "/{ref_organism}/{ref_organism}.done",
#        mt_n = gmap_db_dir + "{ref_organism}/{ref_organism}.done"
    shell:
        "touch {output.mt}"

