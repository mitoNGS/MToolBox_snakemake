#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import click
import collections
import copy
import vcf


def to_list(obj):
    """ If there is only one variant allele, the vcf module will parse it as
        string, not list. This function converts these instances to lists.
    """
    if not isinstance(obj, list):
        obj = [obj]
    return obj


# TODO: this will become one of mtoolbox's subcommands
# TODO: this is ugly, both input and output should be arguments or options
@click.command()
@click.argument("vcf_input")
@click.option("--vcf_output", "-o", default="filtered.vcf",
              help="""Output VCF file.""")
@click.option("--hf_lower", "-l", default=0.20, type=float,
              help="""Lower HF threshold, must be <= 1.0 (Default: 0.20)""")
@click.option("--hf_upper", "-u", default=1.0, type=float,
              help="""Upper HF threshold (Default: 1.0)""")
def HF_filter(vcf_input, vcf_output, hf_lower, hf_upper):
    """ Filter a VCF file based on heteroplasmy frequency. """
    # TODO: check if hf_lower AND HF_UPPER are greater than 1.0 or lower than 0.0
    click.echo("VCF input file: {}".format(vcf_input))
    click.echo("VCF output file: {}".format(vcf_output))
    click.echo("HF lower threshold: {}".format(hf_lower))
    click.echo("HF upper threshold: {}".format(hf_upper))

    vcf_in = vcf.Reader(vcf_input)
    vcf_out = vcf.Writer(open(vcf_output, "w"), vcf_in)

    f_keys = vcf_in.formats.keys()

    for record in vcf_in:
        record_out = copy.deepcopy(record)
        for sx in range(len(record.samples)):
            record_out.samples[sx].data = collections.namedtuple('CallData',
                                                                 f_keys)

            f_vals = [record.samples[sx].data[vx]
                      for vx in range(len(f_keys))]
            data_dict = dict(zip(f_keys, f_vals))
            DP = data_dict['DP']

            try:
                GT = data_dict['GT'].split("/")
                HF_filtered = []
                GT_filtered = []
                CILOW_filtered = []
                CIUP_filtered = []
                SDP_filtered = []
                if sum(to_list(data_dict['HF'])) < (1 - hf_lower):
                    GT_filtered.append(0)
                for x, i in enumerate(to_list(data_dict['HF'])):
                    if i >= hf_lower and i <= hf_upper:
                        HF_filtered.append(i)
                        # if the sample at this POS is heteroplasmic,
                        # len(GT) > 1 and '0' in GT.
                        # GT has REF too so the numbering is shifted
                        if len(GT) > 1:
                            GT_filtered.append(GT[x + 1])
                        # else, grab the GT element with same index
                        else:
                            # could be GT_filtered.append(GT[0]) as well
                            GT_filtered.append(GT[x])
                        CILOW_filtered.append(to_list(data_dict['CILOW'])[x])
                        CIUP_filtered.append(to_list(data_dict['CIUP'])[x])
                        SDP_filtered.append(to_list(data_dict['SDP'])[x])
            except:  # fields might be empty
                HF_filtered = ["."]
                GT_filtered = [".", "."]
                CILOW_filtered = ["."]
                CIUP_filtered = ["."]
                DP = "."
                SDP_filtered = ["."]
            new_vals = ['/'.join([str(n) for n in GT_filtered]),
                        DP, HF_filtered, CILOW_filtered, CIUP_filtered,
                        SDP_filtered]
            record_out.samples[sx].data = record_out.samples[sx].data._make(new_vals)
        # write record only if at least one sample has something left
        ### prints left for testing
        # for sample in vcf_record_new.samples:
        #    print(vcf_record_new.POS, sample, sample.data.GT)
        # if vcf_record_new.POS == 24:
        #    print(set([vcf_record_new.samples[sx].data.GT
        #    for sx in range(len(vcf_record.samples))]))
        ###
        if set([record_out.samples[sx].data.GT
                for sx in range(len(record.samples))]) != {"0"} and \
            set([record_out.samples[sx].data.GT
                 for sx in range(len(record.samples))]) != {"0", "./."}:
            vcf_out.write_record(record_out)

    vcf_out.close()
