# table parsers
def get_other_fields(df, ref_genome_mt, field):
    return list(set(df.loc[df['ref_genome_mt'] == ref_genome_mt, field]))

def get_vcf_files(df, res_dir="results"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{results}/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf".format(results = res_dir, \
                                                                                                                    sample = getattr(row, "sample"), \
                                                                                                                    ref_genome_mt = getattr(row, "ref_genome_mt"), \
                                                                                                                    ref_genome_n = getattr(row, "ref_genome_n")))
    return outpaths

def get_out_files(df, res_dir="results", map_dir="map"):
    outpaths = []
    for row in df.itertuples():
        outpaths.append("{}/OUT_{}_{}_{}/{}/OUT.sam".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n"), map_dir))
    return outpaths

def get_single_vcf_files(df, ref_genome_mt = None):
    ref_genome_mt = ref_genome_mt
    outpaths = []
    for row in df.itertuples():
        if getattr(row, "ref_genome_mt") == ref_genome_mt:
        # "results/OUT_{sample}_{ref_genome_mt}_{ref_genome_n}/vcf.vcf"
            outpaths.append("{}/OUT_{}_{}_{}/vcf.vcf".format(res_dir, getattr(row, "sample"), getattr(row, "ref_genome_mt"), getattr(row, "ref_genome_n")))
    return outpaths

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

############################

def fastqc_raw_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_raw", ext=".fastq.gz"):
    # analysis_tab is a pandas df with a column "sample"
    fastqc_out = []
    for s in analysis_tab["sample"]:
        # fastqc_html_1 and fastqc_html_2 could be strings, but
        # keep them as lists in case one sample has multiple datasets
        fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
        fastqc_out.extend(fastq_files_1)
        fastqc_html_2 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "2", input_folder = infolder)]
        fastqc_out.extend(fastq_files_2)
    return fastqc_out

def fastqc_filtered_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_filtered", ext=".fastq.gz"):
    fastqc_out = []
    for s in analysis_tab["sample"]:
        # fastqc_html_1 and fastqc_html_2 could be strings, but
        # keep them as lists in case one sample has multiple datasets
        #L = set(os.path.join(infolder))
        #glob.glob("{in}/{}*R1{}.format()")
        # expand("results/fastqc_filtered/{sample}")
        # fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
        # fastqc_out.extend(fastq_html_1)
        #######
        for read_type in ["R1", "R2", "U"]:
            fastqc_out.append(os.path.join(outfolder, "{sample}.{read_type}_fastqc.html".format(sample=s, read_type=read_type)))
    return fastqc_out

def fastqc_outputs(analysis_tab = analysis_tab, infolder="data/reads", outfolder="results/fastqc_raw", ext=".fastq.gz", read_types = ["1", "2"]):
    # keyword default values are for raw reads
    fastqc_out = []
    for s in analysis_tab["sample"]:
        # fastqc_html_1 and fastqc_html_2 could be strings, but
        # keep them as lists in case one sample has multiple datasets
        for read_type in read_types:
            for i in read_datasets_inputs(sample = s, read_type = read_type, input_folder = infolder):
                print(i)
                print(os.path.split(i)[1])
                fastqc_html = [os.path.join(outfolder, os.path.split(i)[1].replace(ext, "_fastqc.html"))]
            fastqc_out.extend(fastqc_html)
        # fastqc_html_1 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "1", input_folder = infolder)]
        # fastqc_out.extend(fastq_html_1)
        # fastqc_html_2 = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "2", input_folder = infolder)]
        # fastqc_out.extend(fastq_html_2)
        # fastqc_html_U = [os.path.join(outfolder, i.replace(ext, "_fastqc.html")) for i in read_datasets_inputs(sample = s, read_type = "U", input_folder = infolder)]
        # fastqc_out.extend(fastq_html_U)
    print("fastqc_outputs: {}".format(fastqc_out))
    return fastqc_out

def get_trimmomatic_adapters_path(s):
    trimmomatic_exec_path = s
    return trimmomatic_exec_path.replace("bin/trimmomatic", "share/trimmomatic/adapters/TruSeq3-PE.fa") + ":2:30:10"

