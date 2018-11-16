import pandas as pd
import os

# fields: sample  ref_genome_mt   ref_genome_n
analysis_tab = pd.read_table("analysis.tab", sep = "\t", comment='#')

configfile: "config.yaml"
res_dir = config["results"]
map_dir = config["map_dir"]
log_dir = config["log_dir"]
gsnap_db_dir = config["map"]["gsnap_db_dir"]
# res_dir = "results"
# map_dir = "map"

def get_out_files(df, res_dir="results", map_dir="map"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{}/OUT_{}_{}_{}/{}/OUT.sam".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n"), map_dir))
    return outpaths

def filter_alignments(outmt, outS, outP, OUT, gsnap_db = None):
    sig=1
    pai=1
    print('Reading Results...')
    if sig:
        hgoutsam = outS
        dicsingle={}
        f = open(hgoutsam, 'r')
        for i in f:
            if i.strip()=='': continue
            l=(i.strip()).split('\t')
            if l[2]=='*': continue # the read is not mapped
            # keeping multiple mappings
            if l[0] in dicsingle:
                dicsingle[l[0]].append(l)
            else:
                dicsingle[l[0]]=[l]
        f.close()
    if pai:
        hgoutsam2 = outP
        dicpair={}
        f=open(hgoutsam2)
        for i in f:
            if i.strip()=='': continue
            l=(i.strip()).split('\t')
            if l[2]=='*': continue
            if l[0] in dicpair:
                dicpair[l[0]].append(l)
            else:
                dicpair[l[0]]=[l]
        f.close()

    print('Extracting FASTQ from SAM...')
    mtoutsam = outmt
    dics={}
    f=open(mtoutsam)
    for i in f:
        if i.strip()=='' or i.startswith('@'): continue
        l=(i.strip()).split('\t')
        if l[2]=='*': continue
        if l[0] in dics: dics[l[0]].append(l)
        else: dics[l[0]]=[l]
    f.close()
    
    finalsam = OUT
    out=open(finalsam,'w')
    out.write("@SQ	SN:%s	LN:16569\n" % gsnap_db)
    out.write("@RG	ID:sample	PL:sample	PU:sample	LB:sample	SM:sample\n")
    
    print('Filtering reads...')
    for i in dics:
        ll=dics[i]
        if len(ll)==1: # if the read has one mapping I assume it's SE
            if i in dicsingle:
                r=dicsingle[i] # i is a list of lists (splitted sam lines)
                if len(r)==1: 
                    # check if read aligned on MT when aligned against nuclear+MT
                    # fields checked: RNAME, POS
                    if r[0][2]==ll[0][2] and ll[0][3]==r[0][3]:
                        #good.append('\t'.join(ll[0])+'\n')
                        out.write('\t'.join(ll[0])+'\n')
            else:
                out.write('\t'.join(ll[0])+'\n')
        else:
            if i in dicpair:
                r=dicpair[i]
                if len(r) == 2:
                    if r[0][2]==ll[0][2] and ll[0][3]==r[0][3] and r[1][2]==ll[1][2] and ll[1][3]==r[1][3]:
                        out.write('\t'.join(ll[0])+'\n')
                        out.write('\t'.join(ll[1])+'\n')
            else:
                out.write('\t'.join(ll[0])+'\n')
                out.write('\t'.join(ll[1])+'\n')
    out.close()
    
    print('Outfile saved on %s.' %(finalsam))
    print('Done.')

def gsnap_inputs(wildcards):
    if (seq_type == "pe"):
        return expand("data/reads/{sample}.{strand}.fastq.gz", strand=["R1", "R2"], sample=wildcards.sample)#, filtered_reads=config['proj_dirs']['filtered_reads'])
    elif (seq_type == "se"):
        return expand("data/reads/{sample}.fastq.gz", sample=wildcards.sample)#, filtered_reads=config['proj_dirs']['filtered_reads'])
    elif (seq_type == "both"):
        return expand("data/reads/{sample}{strand}.fastq.gz", strand=[".R1", ".R2", ""], sample=wildcards.sample)#, filtered_reads=config['proj_dirs']['filtered_reads'])

# def gsnap_inputs(wildcards):
#     if (seq_type == "pe"):
#         return expand("{filtered_reads}/{sample}.{strand}.fastq.gz", strand=["R1", "R2"], sample=wildcards.sample, filtered_reads=config['proj_dirs']['filtered_reads'])
#     elif (seq_type == "se"):
#         return expand("{filtered_reads}/{sample}.fastq.gz", sample=wildcards.sample, filtered_reads=config['proj_dirs']['filtered_reads'])
#     elif (seq_type == "both"):
#         return expand("{filtered_reads}/{sample}{strand}.fastq.gz", strand=[".R1", ".R2", ""], sample=wildcards.sample, filtered_reads=config['proj_dirs']['filtered_reads'])

# def filtering_reads_input():
#     if (seq_type == "pe"):
#         return "outhumanP.sam"
#         #return expand("{reads}_{strand}.fastq", strand=["R1", "R2"], reads=wildcards.reads)
#     elif (seq_type == "se"):
#         return "outhumanS.sam"
#         #return expand("{reads}.fastq", reads=wildcards.reads)
#     elif (seq_type == "both"):
#         return "outhumanP.sam", "outhumanS.sam"
#         #return expand("{reads}{strand}.fastq", strand=["_R1", "_R2", ""], reads=wildcards.reads)

seq_type = "both"

outpaths = get_out_files(analysis_tab, res_dir = "results", map_dir = "map")

target_inputs = [
    outpaths ]

rule all:
    input: target_inputs

rule map_MT_PE_SE:
    input:
        #gsnap_inputs,
        R1 = "data/reads/{sample}.R1.fastq.gz",
        R2 = "data/reads/{sample}.R2.fastq.gz",
        SE = "data/reads/{sample}.fastq.gz",
        gsnap_db = gsnap_db_dir + "/{ref_genome_mt}/{ref_genome_mt}.ref081locoffsets64strm"
        #index=gsnap_index
    output:
        outmt_sam = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/outmt.sam"
    params:
        gsnap_db_dir = config["map"]["gsnap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gsnap_db = lambda wildcards, input: os.path.split(input.gsnap_db)[1].split(".")[0],
        # gsnap_db_folder = config['map_exome']['gsnap_db_folder'],
        # gsnap_db = config['map_exome']['gsnap_mt_db'],
        RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample'
    log:
        log_dir + "/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/logmt.txt"
    threads:
        config["map"]["gsnap_threads"]
    run:
        if seq_type == "pe":
            print("PE mode")
            shell("gsnap -D {params.gsnap_db_folder} -d {params.gsnap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} > {output.outmt_sam} 2> {log}")
        if seq_type == "se":
            print("SE mode")
            shell("gsnap -D {params.gsnap_db_folder} -d {params.gsnap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} > {output.outmt_sam} 2> {log}")
        elif seq_type == "both":
            print("PE + SE mode")
            shell("gsnap -D {params.gsnap_db_folder} -d {params.gsnap_db} -A sam --gunzip --nofails --pairmax-dna=500 --query-unk-mismatch=1 {params.RG_tag} -n 1 -Q -O -t {threads} {input[0]} {input[1]} {input[2]} > {output.outmt_sam} 2> {log}")

rule sam2fastq:
    input:
        outmt_sam = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt.sam"
    output:
        outmt1 = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt1.fastq",
        outmt2 = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt2.fastq",
        outmt = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt.fastq"
    # version:
    #     subprocess.getoutput(
    #         "picard SamToFastq --version"
    #         )
    # message:
    #     "Converting SAM files to FASTQ with PicardTools v{version}"
    shell:
        """
        picard SamToFastq \
            I={input.outmt_sam} \
            FASTQ={output.outmt1} \
            SECOND_END_FASTQ={output.outmt2} \
            UNPAIRED_FASTQ={output.outmt}
        """

rule map_nuclear_MT_SE:
    input:
        outmt = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt.fastq",
        gsnap_db = gsnap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.ref081locoffsets64strm"
    output:
        outS = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outS.sam"
    params:
        gsnap_db_dir = config["map"]["gsnap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gsnap_db = lambda wildcards, input: os.path.split(input.gsnap_db)[1].split(".")[0]
        #gsnap_db = config['map']['gsnap_n_mt_db']
    threads:
        config["map"]["gsnap_threads"]
    # version:
    #     subprocess.getoutput(
    #       gsnap --version
    #       )
    #     """
    log:
        logS = log_dir + "/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/logS.sam"
    message:
        "Mapping onto complete human genome (nuclear + mt)... SE reads"
    shell:
        """
        touch {output}
        #gsnap -D {params.gsnap_db_dir} -d {params.gsnap_db} -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt} > {output.outS} 2> {log.logS}
        """

rule map_nuclear_MT_PE:
    input:
        outmt1 = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt1.fastq",
        outmt2 = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt2.fastq",
        gsnap_db = gsnap_db_dir + "/{ref_genome_mt}_{ref_genome_n}/{ref_genome_mt}_{ref_genome_n}.ref081locoffsets64strm"
    output:
        outP = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outP.sam"
    params:
        gsnap_db_dir = config["map"]["gsnap_db_dir"],
        #gsnap_db_folder = config['map']['gsnap_db_folder'],
        gsnap_db = lambda wildcards, input: os.path.split(input.gsnap_db)[1].split(".")[0]
        #gsnap_db = config['map']['gsnap_n_mt_db']
    threads:
        config["map"]["gsnap_threads"]
    # version:
    #     subprocess.getoutput(
    #       gsnap --version
    #       )
    #     """
    log:
        logP = log_dir + "/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/logP.sam"
    message:
        "Mapping onto complete human genome (nuclear + mt)... PE reads"
    shell:
        """
        touch {output}
        #gsnap -D {params.gsnap_db_dir} -d {params.gsnap_db} -A sam --nofails --query-unk-mismatch=1 -O -t {threads} {input.outmt1} {input.outmt2} > {output.outP} 2> {log.logP}
        """

rule filtering_mt_alignments:
    input:
        outmt = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outmt.sam",
        outS = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outS.sam",
        outP = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/outP.sam"
    output:
        sam = "{res_dir}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{map_dir}/OUT.sam"
        #sam = "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/map/OUT.sam"
    params:
        outdir = lambda wildcards, output: os.path.split(output.sam)[0]
    run:
        filter_alignments(input.outmt, \
                          input.outhumanS, \
                          input.outhumanP, \
                          output.sam, \
                          gsnap_db = {params.gsnap_db})
        # try:
        #     os.makedirs(params.outdir)
        # except FileExistsError:
        #     pass
        # with open(output.sam, 'w') as t:
        #     t.write("bwe\n")
