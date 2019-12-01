#!/usr/bin/env python
import click
import csv
import getopt
import os
import os.path
import re
import shutil
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from modules.classifier import tree, NGclassify, consts, parse_mhcs, datatypes
from modules.io_modules import csv, old_table
from modules.bioinf.seqs import SeqList


# folder where to find data for haplogroup classification and functional annotation
# data_file = os.path.dirname(sys.argv[0])

# TODO: what about the data_file argument?

# TODO: this will become one of mtoolbox's subcommands
@click.command()
@click.option("--input_file", "-i", default="mtDNAassembly-contigs.fasta",
              help="""Contig file (Default: mtDNAassembly-contigs.fasta)""")
@click.option("--muscle_exec", "-m", default=shutil.which("muscle"),
              help="""MUSCLE executable path (Default: MUSCLE installation 
              PATH)""")
@click.option("--basename", "-b", default="mtDNAassembly-contigs",
              help="""Basename for output files (Default: mtDNAassembly-contigs)""")
@click.option("--best_results", "-s", default="mt_classification_best_results.csv",
              help="""File with most reliable haplogroup prediction 
              (Default: mt_classification_best_results.csv)""")
def mt_classifier(input_file, muscle_exec, basename, best_results):
    """ Assign haplogroups to contigs and perform functional annotation. """
    pass

# TODO: this seems not working properly, leave it for later


def usage():
    print("""\nAssigns haplogroup to contigs and performs functional annotation
        Options:
        -i		Contig file [mtDNAassembly-Contigs.fasta]
        -m		MUSCLE executable PATH [/usr/local/bin/muscle]
        -b		basename for output files
        -s		file with most reliable haplogroup prediction
        """)


def pickle_csv(csvfile, pickle_fname=None):
    tree_file = csv.reader(open(csvfile, 'rb'))
    if pickle_fname is None:
        pickle_fname = csvfile + '.pickle'
    aplo_list = csv.parse_csv(tree_file)
    htree = tree.HaplogroupTree(aplo_list=aplo_list)
    pickle_file = open(pickle_fname, 'wb')
    pickle_file.write(htree.serialize())


def write_old_table(pickle_fname, out_fname):
    htree = tree.HaplogroupTree(pickle_data=open(pickle_fname, 'rb').read())
    fh = csv.writer(open(out_fname, 'wb'))
    for haplo_name in htree:
        old_table.write_haplogroup(fh, '', htree[haplo_name])


def parse_gmapf9_line(line):
    parsed = line.split('\t')
    last_field = re.findall(r"[\w']+", parsed[2])
    seq_nuc = parsed[1].partition(' ')[2]
    seq_index = parsed[1].partition(' ')[0]
    ref_pos = int(last_field[1])
    ref_nuc = parsed[2][-1]
    return ref_pos, ref_nuc, seq_nuc, seq_index


def parse_gmapf9_file(inhandle):
    contigs_mappings = [[]]
    h = inhandle.readlines()
    c = 0
    mutations = []
    while c < len(h):
        # end coordinate of last contig
        if c == len(h)-1:
            contigs_mappings[-1].append(parse_gmapf9_line(h[c])[0])
        if h[c][0] != '>':
            ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
            # insertion
            if ref_nuc == ' ' and seq_nuc != ' ':
                # gmap assigns the position of the next nucleotide to the insertion
                pos_ins = ref_pos - 1
                ins = [seq_nuc]
                c += 1
                ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
                while c < len(h) and (ref_nuc == ' ' and seq_nuc != ' '):
                    ins.append(seq_nuc)
                    c += 1
                    ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
                mut = datatypes.Insertion("%d.%s" % (pos_ins, ''.join(ins)))
                mutations.append(mut)
            # deletion
            elif ref_nuc != ' ' and seq_nuc == ' ':
                pos_del = ref_pos
                c += 1
                ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
                while c < len(h) and (ref_nuc != ' ' and seq_nuc == ' '):
                    c += 1
                    ref_pos, ref_nuc, seq_nuc, seq_index = parse_gmapf9_line(h[c])
                if pos_del == ref_pos-1:
                    print("{}d".format(pos_del))
                    mut = datatypes.Deletion("%dd" % pos_del)
                    mutations.append(mut)
                else:
                    print("{}-{}d".format(pos_del, ref_pos-1))
                    mut = datatypes.Deletion("%d-%dd" % (pos_del, ref_pos-1))
                    mutations.append(mut)
            # mismatch
            elif ref_nuc != seq_nuc:
                if seq_nuc != 'N':
                    # Transition
                    if ((ref_nuc in consts.PUR and seq_nuc in consts.PUR) or
                            (ref_nuc in consts.PYR and seq_nuc in consts.PYR)):
                        print("{}{}".format(ref_pos, seq_nuc))
                        mut = datatypes.Transition(ref_pos)
                        mutations.append(mut)
                    # Transversion
                    if ((ref_nuc in consts.PUR and seq_nuc in consts.PYR) or
                            (ref_nuc in consts.PYR and seq_nuc in consts.PUR)):
                        mut = datatypes.Transversion("%d%s" % (ref_pos, seq_nuc))
                        mutations.append(mut)
                c += 1
            else:
                c += 1
        else:
            # first contig
            if len(contigs_mappings) == 1 and len(contigs_mappings[-1]) == 0:
                contigs_mappings[-1].append(parse_gmapf9_line(h[c+1])[0])
            # all the others
            else:
                contigs_mappings[-1].append(parse_gmapf9_line(h[c-1])[0])
                contigs_mappings.append([parse_gmapf9_line(h[c+1])[0]])
            c += 1
    # don't know if contig coordinate sorting is needed but I'll do anyway
    contigs_mappings.sort()
    return mutations, contigs_mappings


def merge_tables(f, g, h):
    fgh = f + g + h
    mergedlist = []
    for jj in fgh:
        if jj not in mergedlist:
            mergedlist.append(jj)
    o = []
    o.append(["", "RSRS", "MHCS", "rCRS"])
    y = "yes"
    n = ""
    for i in mergedlist:
        if i in f and i in g and i in h:
            o.append([i.pprint(), y, y, y])
        elif i in f and i in g:
            o.append([i.pprint(), y, y, n])
        elif i in f and i in h:
            o.append([i.pprint(), y, n, y])
        elif i in g and i in h:
            o.append([i.pprint(), n, y, y])
        elif i in f:
            o.append([i.pprint(), y, n, n])
        elif i in g:
            o.append([i.pprint(), n, y, n])
        elif i in h:
            o.append([i.pprint(), n, n, y])
    return o


def align_sequence(muscle_exe, obj=None, rif=None):
    """obj is a Bio.SeqRecord.SeqRecord"""
    if rif is None:
        rif = SeqRecord(Seq(consts.RCRS), id='RSRS', name='RSRS')
    seq_diff = NGclassify.SequenceDiff()
    seq_diff.gen_diff(muscle_exe=muscle_exe, rif=rif, obj=obj)
    return seq_diff


def h_analysis(htrees, seq_diff, regions, mhcs_dict):
    a = NGclassify.Classify()
    for htree, name in htrees:
        print("Classification according to tree: {}".format(name))
        a.classify_by_tree(htree, seq_diff, regions)
        print("genome_state is {}".format(a.get_genome_state()))
        (haplo_stats_sorted, haplo_best) = a.prediction_sorting()
        print(haplo_best)
        print("="*20)
        # TODO: mhcss is not used anywhere
        mhcss = a.get_mhcss(mhcs_dict)
    print('-'*30)
    return a


def load_sequences(fname):
    a = SeqList()
    a.load_file(fname)
    print("Loaded {} contig sequences".format(len(a)))
    return a


# TODO: seq_diff, seq_diff_mhcs, seq_diff_rcrs are not used anywhere
def write_output(class_obj, seq_diff, seq_diff_mhcs, seq_diff_rcrs, merged_tables, outfile):
    print("Writing results for sequence {}".format(outfile))
    class_obj.pprint(open(outfile + '.csv', 'w'))
    class_obj.pprint_sorted(open(outfile + '.sorted.csv', 'w'))
    merged_tables_file = open(outfile + '_merged_diff.csv', 'w')
    for row in merged_tables:
        merged_tables_file.write(','.join(row)+'\n')


def main_mt_hpred(contig_file='mtDNAassembly-contigs.fasta',
                  muscle_exe="/usr/bin/muscle",
                  basename="mtDNAassembly-contigs",
                  best_results_file='mt_classification_best_results.csv',
                  data_file=None):
    print("Your best results file is {}".format(best_results_file))
    # sample name
    f = os.path.abspath(contig_file)
    sample_name = contig_file.split('-')[0]
    # haplogroup tree parsing
    htrees = [
        (tree.HaplogroupTree(
            pickle_data=open(os.path.join(data_file, 'phylotree_r17.pickle'),
                             'rb').read()
        ),
         os.path.join(data_file, 'phylotree_r17.pickle'))
    ]
    print(htrees)
    # mhcs parsing
    mhcs_dict = parse_mhcs.parse2mhcs_dict(os.path.join(data_file, 'mhcs.tab'))
    
    print("\nLoading contig sequences from file {}".format(contig_file))
    contig_array = SeqIO.index(contig_file, 'fasta')

    print("\nAligning Contigs to mtDNA reference genome...\n")
    
    # update each contig's SeqDiff
    for x, contig in enumerate(list(contig_array.keys())):
        if x == 0:
            contig_seq_diff = align_sequence(muscle_exe,
                                             obj=contig_array[contig],
                                             rif=SeqRecord(Seq(consts.RCRS),
                                                           id='RSRS',
                                                           name='RSRS'))
            # avoid having long gaps at 5' and 3' (not actual gaps but due to the alignment)
            contig_seq_diff.find_segment()
            contig_seq_diff.regions.append([contig_seq_diff.start,
                                            contig_seq_diff.end])
        else:
            incoming_seqdiff = align_sequence(muscle_exe,
                                              obj=contig_array[contig],
                                              rif=SeqRecord(Seq(consts.RCRS),
                                                            id='RSRS',
                                                            name='RSRS'))
            incoming_seqdiff.find_segment()
            contig_seq_diff.diff_list.extend(incoming_seqdiff.diff_list)
            contig_seq_diff.regions.append([incoming_seqdiff.start,
                                            incoming_seqdiff.end])
    
    print("\nSequence haplogroup assignment\n")
    seq_classify = h_analysis(htrees, contig_seq_diff,
                              contig_seq_diff.regions, mhcs_dict)
    seq_classify.sample_name = sample_name
    
    print("Contig alignment to MHCS and rCRS")
    m = list(seq_classify.mhcss)[0]
    print("Aligning contigs to MHCS SeqDiff object")
    its_mhcs = SeqRecord(Seq(mhcs_dict[m]), id = m, name = m)
    for x, contig in enumerate(contig_array):
        if x == 0:
            contig_mhcs_seq_diff = align_sequence(muscle_exe,
                                                  obj=contig_array[contig],
                                                  rif=its_mhcs)
            contig_mhcs_seq_diff.find_segment()
            contig_mhcs_seq_diff.regions.append([contig_seq_diff.start,
                                                 contig_seq_diff.end])
        else:
            incoming_mhcs_seqdiff = align_sequence(muscle_exe,
                                                   obj=contig_array[contig],
                                                   rif=its_mhcs)
            incoming_mhcs_seqdiff.find_segment()
            contig_mhcs_seq_diff.diff_list.extend(incoming_mhcs_seqdiff.diff_list)
            contig_mhcs_seq_diff.regions.append([incoming_mhcs_seqdiff.start,
                                                 incoming_mhcs_seqdiff.end])
    
    print("rCRS SeqDiff object")
    # TODO: rcrs here is useless
    rcrs = datatypes.Sequence('rCRS', consts.rcrs)
    rcrs = SeqRecord(Seq(consts.rcrs), id=m, name=m)
    for x, contig in enumerate(contig_array):
        if x == 0:
            contig_rcrs_seq_diff = align_sequence(muscle_exe,
                                                  obj=contig_array[contig],
                                                  rif=rcrs)
            contig_rcrs_seq_diff.find_segment()
            contig_rcrs_seq_diff.regions.append([contig_seq_diff.start,
                                                 contig_seq_diff.end])
        else:
            incoming_rcrs_seqdiff = align_sequence(muscle_exe,
                                                   obj=contig_array[contig],
                                                   rif=rcrs)
            incoming_rcrs_seqdiff.find_segment()
            contig_rcrs_seq_diff.diff_list.extend(incoming_rcrs_seqdiff.diff_list)
            contig_rcrs_seq_diff.regions.append([incoming_rcrs_seqdiff.start,
                                                 incoming_rcrs_seqdiff.end])
    
    # try gathering diff from reference sequences
    print("Merging seq_diffs...")
    mergedtables = merge_tables(contig_seq_diff.diff_list,
                                contig_mhcs_seq_diff.diff_list,
                                contig_rcrs_seq_diff.diff_list)
    # OUTPUTS
    open(best_results_file, 'a').write(','.join(
        [basename, ';'.join([i[0] for i in seq_classify.haplo_best.items()])]
    ) + '\n')
    return (seq_classify, contig_seq_diff, contig_mhcs_seq_diff,
            contig_rcrs_seq_diff, mergedtables)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:m:b:s:d:")
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit()
    contig_file = 'mtDNAassembly-contigs.fasta'
    muscle_exe = shutil.which('muscle')
    basename = 'mtDNAassembly-contigs'
    best_results_file = 'mt_classification_best_results.csv'
    for o, a in opts:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-d":
            data_file = a
        elif o == "-i":
            contig_file = a
        elif o == "-m":
            muscle_exe = a
        elif o == "-b":
            basename = a
        elif o == "-s":
            best_results_file = a
        else:
            assert False, "Unhandled option."

    # TODO: what is data_file?
    if data_file is None:
        sys.exit(("You must specify the folder where data for mt_classifier "
                  "execution are located. Abort."))
    else:
        # TODO: asserting like this is not good
        assert os.path.exists(os.path.join(data_file, "phylotree_r17.pickle"))
    (sc, contig_seq_diff, contig_mhcs_seq_diff,
     contig_rcrs_seq_diff, mergedtables) = main_mt_hpred(contig_file=contig_file,
                                                         muscle_exe=muscle_exe,
                                                         basename=basename,
                                                         best_results_file=best_results_file,
                                                         data_file=data_file)
    # TODO: what is seq_classify? guess it may be "sc" defined above..?
    write_output(seq_classify, contig_seq_diff.diff_list,
                 contig_mhcs_seq_diff.diff_list, contig_rcrs_seq_diff.diff_list,
                 mergedtables, basename)
