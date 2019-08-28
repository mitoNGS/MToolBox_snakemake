#!/usr/bin/env python

import collections, copy, vcf, argparse, sys

def to_list(obj):
    """
    If there is only one variant allele, the vcf module will parse it as string,
    not list. This function converts these instances in lists.
    """
    if type(obj) is not list:
        obj = [obj]
    return obj

#hf_threshold = 0.20

parser = argparse.ArgumentParser(description='Filter VCF file based on heteroplasmy frequency.', \
                                    add_help = False, \
                                    formatter_class=argparse.RawTextHelpFormatter)

group_required = parser.add_argument_group('required arguments')
group_required.add_argument('--vcf-input', '-i', action="store", dest="vcf_input", help = "VCF file to be filtered")
group_optional = parser.add_argument_group('optional arguments')
group_optional.add_argument('--vcf-output', '-o', action="store", dest="vcf_output", help = "(default: <vcf_input>.<hf_threshold>.vcf)")
group_optional.add_argument('--hf-threshold', '-t', action="store", dest="hf_threshold", type=float, default=0.20, help = "Must be <= 1 (default: %(default)s)")
group_optional.add_argument("-h", "--help", action="help", help="show this help message and exit")

args = parser.parse_args()

if args.vcf_input == None:
    parser.print_help()
    sys.exit("\nERROR: --vcf-input is required. Exit.\n")
else:
    vcf_input = args.vcf_input

if args.hf_threshold > 1:
    parser.print_help()
    sys.exit("\nERROR: --hf-threshold must be <= 1. Exit.\n")
else:
    hf_threshold = args.hf_threshold
    
if args.vcf_output == None:
    vcf_output = args.vcf_input.replace(".vcf", "_filt{}.vcf".format(hf_threshold))
    print("VCF output is {}".format(vcf_output))

VCF = vcf.Reader(filename = vcf_input)
VCF_new = vcf.Writer(open(vcf_output, 'w'), VCF)

f_keys = VCF.formats.keys()

for vcf_record in VCF:
    vcf_record_new = copy.deepcopy(vcf_record)
    for sx in range(len(vcf_record.samples)):
        vcf_record_new.samples[sx].data = collections.namedtuple('CallData', f_keys)
        #for vx in f_keys:
            #print(vx, getattr(vcf_record.samples[sx].data, vx))
        f_vals = [vcf_record.samples[sx].data[vx] for vx in range(len(f_keys))]
        data_dict = dict(zip(f_keys, f_vals))
        DP = data_dict['DP']
        ##########
        try:
            GT = data_dict['GT'].split("/")
            HF_filtered = []
            GT_filtered = []
            CILOW_filtered = []
            CIUP_filtered = []
            if sum(to_list(data_dict['HF'])) < 1 - hf_threshold:
                GT_filtered.append(0)
            for x, i in enumerate(to_list(data_dict['HF'])):
                if i >= hf_threshold:
                    HF_filtered.append(i)
                    GT_filtered.append(GT[x+1]) # GT has REF too so the numbering is shifted
                    CILOW_filtered.append(to_list(data_dict['CILOW'])[x])
                    CIUP_filtered.append(to_list(data_dict['CIUP'])[x])
        except: # fields might be empty
            HF_filtered = ["."]
            GT_filtered = [".", "."]
            CILOW_filtered = ["."]
            CIUP_filtered = ["."]
            DP = "."
        new_vals = ['/'.join([str(n) for n in GT_filtered]), \
                   DP, \
                   HF_filtered, \
                   CILOW_filtered, \
                   CIUP_filtered]
        vcf_record_new.samples[sx].data = vcf_record_new.samples[sx].data._make(new_vals)
    # write record only if at least one sample has something left
    ### prints left for testing
    #for sample in vcf_record_new.samples:
    #    print(vcf_record_new.POS, sample, sample.data.GT)
    #if vcf_record_new.POS == 24:
    #    print(set([vcf_record_new.samples[sx].data.GT for sx in range(len(vcf_record.samples))]))
    ###
    if set([vcf_record_new.samples[sx].data.GT for sx in range(len(vcf_record.samples))]) != set(["0"]) and \
        set([vcf_record_new.samples[sx].data.GT for sx in range(len(vcf_record.samples))]) != set(["0", "./."]):
        VCF_new.write_record(vcf_record_new)

VCF_new.close()
