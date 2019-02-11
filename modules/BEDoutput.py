import sys, os, glob, math, gzip
#print sys.version
import re
import ast
from collections import OrderedDict
import vcf

def get_BED_score(VCF_entry, x):
    """
    HF expressed as %
    """
    return VCF_entry._sample_indexes[0][1][2][x]*100

def assignRGB(mut_type):
    if mut_type == "ins": # green
        return "0,255,0"
    elif mut_type == "del": # red
        return "255,0,0"
    else: # blue, maybe add a stricter control as "mism"
        return "0,0,255"

def write_bed_header(bed_handle, seq_name = "seq"):
    bed_handle.write('track name="{s} mt variants" description="{s} mitochondrial variants" viewLimits=0:100 useScore=1 itemRgb="On"\n'.format(s=seq_name))


def BEDoutput(VCF_RECORDS, seq_name="seq", bedfile="bed.bed"):
    outBED = open(bedfile, 'w')
    write_bed_header(outBED, seq_name=seq_name)
    for r in VCF_RECORDS:
        print(r.__dict__)
        for x, allele in enumerate(r.ALT):
            if r.ALT[x] == r.REF[0]: continue
            START = r.POS + (len(r.REF[0]) - 1) - 1
            END = START + 1
            SCORE = get_BED_score(r, x)
            if r.TYPEVAR[x] == 'ins':
                # get starting position:
                # - get the one from VCF, but it can be long > 1 since it could be the same for a deletion
                # - add the length of the REF allele-1 to get the pos of the last nt in REF
                # - subtract 1 because START in BED is 0-based
                # get end position:
                # - start position + 1
                NAME = "{}.{}".format(END, allele.replace(r.REF[0], '', 1))
                #SCORE = get_BED_score(r, x)
            elif r.TYPEVAR[x] == 'del':
                #START = r.POS - 1
                #START = r.POS + (len(r.REF[0]) - 1) - 1
                #END = START + 1
                #END = r.POS + (len(r.REF[0]) - 1) - 1
                if START+1 == END:
                    NAME = "{}d".format(END)
                else:
                    NAME = "{}-{}d".format(START+1, END)
                #SCORE = get_BED_score(r, x)
            else:
                NAME = "{}{}".format(r.POS, r.ALT[x])
            outBED.write("\t".join(map(lambda x: str(x), [r.CHROM, \
                                           START, \
                                           END, \
                                           NAME, \
                                           #r.REF[0], \
                                           SCORE, \
                                           ".", \
                                           START, \
                                           END, \
                                           assignRGB(r.TYPEVAR[x])
    #                                       allele, \
    #                                       r.TYPEVAR[x]
                                                  ]))+"\n")
    #         print("\t".join(map(lambda x: str(x), [r.CHROM, \
    #                                        START, \
    #                                        END, \
    #                                        NAME, \
    #                                        #r.REF[0], \
    #                                        SCORE, \
    #                                        ".", \
    #                                        START, \
    #                                        END, \
    #                                        assignRGB(r.TYPEVAR[x])
    # #                                       allele, \
    # #                                       r.TYPEVAR[x]
    #                                               ])))
    outBED.close()

if __name__ == '__main__':
	print("This script is used only when called by MToolBox.")
	pass
