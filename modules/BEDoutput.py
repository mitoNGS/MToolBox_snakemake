import sys, os, glob, math, gzip
#print sys.version
import re
import ast
from collections import OrderedDict
import vcf
import pandas as pd
from modules.mtVariantCaller import get_consensus_single

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
        #print(r.__dict__)
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

def FASTAoutput(vcf_dict = None, ref_mt = None, contigs = [], fasta_out = "fasta.fasta", hf = 0.8):
    mut_events = vcf_dict
    crf = True
    if crf: f=open(fasta_out,'w')
    x=1
    for i in contigs:
        #initialize new_i
        new_i = i
        #write fasta header
        f.write('>Contig.%i|%i-%i\n' %(x,new_i[0][0],new_i[0][1]))
        #print "A contig, ", i
        if crf:
            string_seq = i[1]
            #print "String seq is", string_seq
            nuc_index = i[0][0]
            dict_seq = {}
            # the sequence string at 
            for nuc in string_seq:
                dict_seq[nuc_index] = nuc
                nuc_index += 1
            #print "original dict_seq is", dict_seq
            # add info for consensus dictionary
            consensus_single = get_consensus_single(mut_events[list(mut_events.keys())[0]],hf=hf)
            #print consensus_single
            # alter dict_seq keys for the implementation
            # of the consensus information
            #
            #print "CONSENSUS SINGLE: ", consensus_single
            #check if there are repeated positions with different mut type
            if len(consensus_single) == 0:
                print('no variants found in this contig {0}\n').format(x)
                pass
            else:
                df= pd.DataFrame(consensus_single)
                positions=df[0]
                dup_positions = positions[positions.duplicated()].values
                for x in dup_positions:
                    d = df[df[0]==x][2] #check the mut type. If ins, report ins instead of del or mism 
                    if 'ins' in d.values:
                        idx = d[d!='ins'].index[0]
                        df.drop(df.index[[idx]],inplace=True)
                    elif 'del' in d.values:
                        idx = d[d!='del'].index[0]
                        df.drop(df.index[[idx]],inplace=True) #If ambiguity between mism and del, report deletion instead of mism in the consensus
                for idx in df.index:
                    if df[0][idx] in dict_seq.keys(): #if position is in the dict
                        if df[2][idx] == 'mism': #if mut type is mism
                            dict_seq[df[0][idx]] = df[1][idx][0] #then substitute the dict value with the correspondent nt sequence
                        elif df[2][idx] == 'ins':
                            dict_seq[df[0][idx]] = df[1][idx][0]
                        elif df[2][idx] == 'del':
                            for deleted_pos in df[1][idx]:
                                if deleted_pos in dict_seq:
                                    del(dict_seq[deleted_pos])
                                else:
                                    pass #do not try to delete the position from the contig as the position is not present in it (this happens when the deletion is downstream to the end of the contig!
                # sort positions in dict_seq and join to have the sequence
                contig_seq = ''
                #print "dict_seq is", dict_seq.keys()
                for j in sorted(dict_seq.keys()):
                    contig_seq += dict_seq[j]
                #print contig_seq
                new_i = ((i[0][0], i[0][1]), contig_seq)
                #contigs_wdict.append(new_i)
                #f.write('>Contig.%i|%i-%i\n' %(x,new_i[0][0],new_i[0][1]))
                #f.write('>Contig.%i|%i-%i\n' %(x,i[0][0],i[0][1]))
        #dass[i[0]]=[0,0,0,0,0]
        for j in range(0,len(new_i[1]),60):
            if crf:
                f.write(new_i[1][j:j+60]+'\n')
        x+=1
    if crf: f.close()
    #pass

if __name__ == '__main__':
    print("This script is used only when called by MToolBox.")
    pass
