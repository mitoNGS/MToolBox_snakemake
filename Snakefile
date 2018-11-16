import pandas as pd
import os

#analysis_tab = pd.read_table(config["samples_data"], comment='#').set_index("sample", drop=False)
# fields: sample  ref_genome_mt   ref_genome_n
analysis_tab = pd.read_table("analysis.tab", sep = "\t", comment='#')

def get_out_files(df, res_dir="results", map_dir="map"):
    outpaths = []
    for row in df.itertuples():
        # print(row)
        # print(getattr(row, "sample"))
        # print(getattr(row, "ref_genome_mt"))
        # print(getattr(row, "ref_genome_n"))
        outpaths.append("{}/OUT_{}_{}_{}/{}/OUT.sam".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n"), map_dir))
        #outpaths.append("OUT_{}_{}_{}/{}/OUT.sam".format(row["sample"], row["ref_genome_mt"], row["ref_genome_n"], map_dir))
    return outpaths

outpaths = get_out_files(analysis_tab, res_dir = "results", map_dir = "map")

target_inputs = [
    outpaths ]

rule all:
    input: target_inputs

rule map:
    output:
        sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/OUT.sam"
    # params:
    #     outdir = os.path.split(output.outpaths)[0]
    run:
        """
        outdir=os.path.split(output.sam)[0]
        os.makedirs(outdir)
        t = open(outdir + "/" + output.sam, 'w')
        t.close()
        """