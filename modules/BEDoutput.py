#!/usr/bin/env python
from typing import List, Optional

import pandas as pd

from modules.mtVariantCaller import get_consensus_single


# TODO: change module name


def assign_rgb(mut_type: Optional[str] = None) -> str:
    """ Return the RGB color code related to the given mutation type.

    Parameters
    ----------
    mut_type : str, None
        Type of mutation [ins, del] (default: None).

    Returns
    -------
    str
        RGB color code.
    """
    if mut_type == "ins":  # green
        return "0,255,0"
    elif mut_type == "del":  # red
        return "255,0,0"
    else:  # blue, maybe add a stricter control as "mism"
        return "0,0,255"


# TODO: check docstring
def get_bed_score(entry, x):
    """ Return HF as a percentage.

    Parameters
    ----------
    entry :
        VCF entry to process?
    x :
        Who knows.

    Returns
    -------
    clueless
        Don't know.
    """
    print(entry._sample_indexes)
    return entry._sample_indexes[0][1][2][x] * 100


def create_bed_header(seq_name: str = "seq") -> str:
    """ Create the header for a BED file.

    Parameters
    ----------
    seq_name : str
        Sequence name?

    Returns
    -------
    str
        Formatted BED header.
    """
    header = ('track name="{s} mt variants" '
              'description="{s} mitochondrial variants" '
              'viewLimits=0:100 useScore=1 '
              'itemRgb="On"\n').format(s=seq_name)
    return header


# def write_bed_header(bed_handle, seq_name: str = "seq"):
#     """ Write the header for the given BED file.
#
#     Parameters
#     ----------
#     bed_handle :
#         Input BED file.
#     seq_name : str
#         Sequence name?
#     """
#     bed_handle.write(('track name="{s} mt variants" '
#                       'description="{s} mitochondrial variants" '
#                       'viewLimits=0:100 useScore=1 '
#                       'itemRgb="On"\n').format(s=seq_name))


# TODO: check docstring
def bed_output(vcf_records: List, seq_name: str = "seq",
               bedfile: str = "bed.bed"):
    """ Write a set of VCF records to a BED file.

    Parameters
    ----------
    vcf_records : List
        List of VCF records to write to BED.
    seq_name : str
        Sequence name? (default: seq).
    bedfile : str
        Output BED file name (default: bed.bed).
    """
    with open(bedfile, "w") as out_bed:
        header = create_bed_header(seq_name=seq_name)
        out_bed.write(header)
        for r in vcf_records:
            for x, allele in enumerate(r.ALT):
                if r.ALT[x] == r.REF[0]:
                    continue

                START = r.POS + (len(r.REF[0]) - 1) - 1
                END = START + 1
                print(r)
                #SCORE = get_bed_score(r, x) #TODO I change this because something is wrong with the get_bed_score function. Need to fix that first
                SCORE = 0

                if r.TYPEVAR[x] == 'ins':
                    # TODO: this goes into the docstring
                    # get starting position:
                    # - get the one from VCF, but it can be long > 1
                    #   since it could be the same for a deletion
                    # - add the length of the REF allele-1 to get the pos of
                    #   the last nt in REF
                    # - subtract 1 because START in BED is 0-based
                    # get end position:
                    # - start position + 1
                    NAME = "{}.{}".format(END, allele.replace(r.REF[0], '', 1))

                elif r.TYPEVAR[x] == 'del':
                    if START + 1 == END:
                        NAME = "{}d".format(END)
                    else:
                        NAME = "{}-{}d".format(START + 1, END)

                else:
                    NAME = "{}{}".format(r.POS, r.ALT[x])

                str_y = map(lambda y: str(y),
                            [r.CHROM, START, END, NAME, SCORE, ".", START, END,
                             assign_rgb(r.TYPEVAR[x])])
                out_bed.write("\t".join(str_y) + "\n")


# TODO: ref_mt is not used, is it necessary?
# TODO: better to switch position of hf_max and hf_min
def fasta_output(vcf_dict: Optional[dict] = None, ref_mt=None,
                 contigs: Optional[List] = None,
                 fasta_out: str = "fasta.fasta",
                 hf_max: float = 0.8, hf_min: float = 0.2):
    if contigs is None:
        contigs = []
    mut_events = vcf_dict
    # the following does not make sense as crf is always True,
    # converting to with clause
    # crf = True
    # if crf:
    #     f = open(fasta_out, 'w')
    with open(fasta_out, "w") as f:
        x = 1
        for i in contigs:
            # initialize new_i
            new_i = i
            # write fasta header
            f.write(">Contig.{}|{}-{}\n".format(x, new_i[0][0], new_i[0][1]))

            # if crf:
            string_seq = i[1]
            nuc_index = i[0][0]
            dict_seq = {}
            # the sequence string at
            for nuc in string_seq:
                dict_seq[nuc_index] = nuc
                nuc_index += 1

            # add info for consensus dictionary
            consensus_single = get_consensus_single(
                mut_events[list(mut_events.keys())[0]],
                hf_max=hf_max, hf_min=hf_min
            )

            # alter dict_seq keys for the implementation
            # of the consensus information

            # check if there are repeated positions with different mut type
            if len(consensus_single) == 0:
                print('no variants found in this contig {}\n'.format(x))
                pass
            else:
                df = pd.DataFrame(consensus_single)
                positions = df[0]
                dup_positions = positions[positions.duplicated()].values
                for x in dup_positions:
                    # check the mut type. If ins, report ins instead of del or mism
                    d = df[df[0] == x][2]
                    if 'ins' in d.values:
                        idx = d[d != 'ins'].index[0]
                        df.drop(df.index[[idx]], inplace=True)
                    elif 'del' in d.values:
                        idx = d[d != 'del'].index[0]
                        # If ambiguity between mism and del, report deletion
                        # instead of mism in the consensus
                        df.drop(df.index[[idx]], inplace=True)
                for idx in df.index:
                    if df[0][idx] in dict_seq.keys():  # if position is in the dict
                        if df[2][idx] == 'mism':  # if mut type is mism
                            # then substitute the dict value with the
                            # correspondent nt sequence
                            dict_seq[df[0][idx]] = df[1][idx][0]
                        elif df[2][idx] == 'ins':
                            dict_seq[df[0][idx]] = df[1][idx][0]
                        elif df[2][idx] == 'del':
                            for deleted_pos in df[1][idx]:
                                if deleted_pos in dict_seq:
                                    del(dict_seq[deleted_pos])
                                else:
                                    # do not try to delete the position from the contig as the
                                    # position is not present in it (this happens when the
                                    # deletion is downstream to the end of the contig!)
                                    pass
                # sort positions in dict_seq and join to have the sequence
                contig_seq = ''

                for j in sorted(dict_seq.keys()):
                    contig_seq += dict_seq[j]

                new_i = ((i[0][0], i[0][1]), contig_seq)

            for j in range(0, len(new_i[1]), 60):
                # if crf:
                f.write(new_i[1][j: j+60] + '\n')
            x += 1
    # if crf:
    #     f.close()


if __name__ == '__main__':
    print("This script is used only when called by MToolBox.")
    pass
