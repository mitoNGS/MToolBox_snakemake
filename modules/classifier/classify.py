#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from collections import OrderedDict
import subprocess
import sys

from .consts import ambiguity, DNA, PUR, PYR, RCRS
from .datatypes import (
    Deletion, Insertion, Sequence, Transition, Transversion, Unknown
)
from .parse_mhcs import which_mhcs_lite

rsrs_seq = RCRS  # actually it's RSRS


def get_snps_from_vcf_dict(vcf_dict, ref_seq=rsrs_seq):
    mutations = []
    for snp in vcf_dict:
        if snp[-1] == 'del':
            # [8270, [8271, 8272, 8273, 8274, 8275, 8276, 8277, 8278, 8279], 'del']
            mut = Deletion("%d-%dd" % (snp[1][0], snp[1][-1]))
        elif snp[-1] == 'mism':
            # [73, ['G'], 'mism'] ma in realtà vorrei
            # [73, ['A', 'G'], 'mism'] dove 'A' sarebbe il nt in RSRS, 'G' la mutazione
            ref = ref_seq[snp[0]-1]
            var = snp[1][0]

            if (ref in PUR and var in PUR) or (ref in PYR and var in PYR):
                mut = Transition(snp[0])
            elif (ref in PUR and var in PYR) or (ref in PYR and var in PUR):
                mut = Transversion("%d%c" % (snp[0], var))
            else:
                mut = Unknown(snp[0])
                mut.change = var
        elif snp[-1] == 'ins':
            # [309, ['CCT'], 'ins']
            mut = Insertion("%d.%s" % (snp[0], snp[1]))
            pass
        mutations.append(mut)
    return mutations


def get_snps(rif, inc, start_pos=0, gap='-'):
    # pos_a è quella assoluta, relativa cioè alle 2 sequenze allineate, n_gaps si
    # riferisce al numero di gap presenti fino a quel punto in anderson
    pos_a = n_gaps = start_pos
    alg_len = len(rif)
    mutations = []
    while pos_a < (alg_len + start_pos):
        x = rif[pos_a]
        y = inc[pos_a]
        if x != y:
            if x != gap and y != gap:
                # SNP
                # Transizione
                if (x in PUR and y in PUR) or (x in PYR and y in PYR):
                    mut = Transition(pos_a-n_gaps + 1)
                    # Nel caso il genoma di riferimento non sia Anderson
                    mut.change = y
                    mutations.append(mut)
                # Trasversione
                elif (x in PUR and y in PYR) or (x in PYR and y in PUR):
                    mut = Transversion("%d%c" % (pos_a-n_gaps + 1, y))
                    mutations.append(mut)
                # Ambiguity
                elif y in ambiguity.keys():
                    for i in ambiguity[y]:
                        # retain mutation defined by ambiguity only if it's
                        # not equal to ref sequence
                        if i != x:
                            if (x in PUR and i in PUR) or (x in PYR and i in PYR):
                                mut = Transition(pos_a - n_gaps + 1)
                                # Nel caso il genoma di riferimento non sia Anderson
                                mut.change = i
                                mut.ambiguity = y
                                mutations.append(mut)
                            elif (x in PUR and i in PYR) or (x in PYR and i in PUR):
                                mut = Transversion("%d%c" % (pos_a-n_gaps + 1, i))
                                mut.ambiguity = y
                                mutations.append(mut)
                # Non identificabile: N o altre ambiguità
                else:
                    mut = Unknown(pos_a-n_gaps+1)
                    mut.change = y
                    mutations.append(mut)
                
                pos_a += 1
            elif y != gap and x == gap:
                # Inserzione
                pos_i = pos_a - n_gaps
                ins_seq = [y]
                pos_a += 1
                n_gaps += 1
                try:
                    x = rif[pos_a]
                    y = inc[pos_a]
                except IndexError:
                    # caso limite: l'inserzione e' di lunghezza 1 alla fine
                    # dell'allineamento
                    x = rif[pos_a-1]
                    y = inc[pos_a-1]
                while pos_a < alg_len-1 and ((x == gap and y != gap) or (x == y == gap)):
                    if y != gap:
                        ins_seq.append(y)
                    pos_a += 1
                    n_gaps += 1
                    x = rif[pos_a]
                    y = inc[pos_a]
                if pos_a == alg_len - 1:
                    pos_a += 1
                mut = Insertion("%d.%s" % (pos_i, ''.join(ins_seq)))
                mutations.append(mut)
            elif y == gap and x != gap:
                # Delezione
                pos_d = pos_a-n_gaps+1
                pos_a += 1
                if pos_a < alg_len:
                    x = rif[pos_a]
                    y = inc[pos_a]
                    while pos_a < alg_len-1 and ((y == gap and x != gap) or (x == y == gap)):
                        pos_a += 1
                        if x == y == gap: n_gaps += 1
                        x = rif[pos_a]
                        y = inc[pos_a]
                    if pos_a == alg_len - 1: pos_a += 1
                mut = Deletion("%d-%dd" % (pos_d, pos_a-n_gaps))
                mutations.append(mut)
        else:
            # accrocchio per permettere di associare quelle riconosciute come retromutazioni nel'albero con quelle che poi in rCRS sono presenti
            pos_a += 1
            # basta controllarne uno, si sa che sono uguali
            if x == gap:
                n_gaps += 1
    return mutations


def compare_mutations(h_pos_list, s_pos_list, start=None, end=None):
    # h_pos_list is the list of variants defining a haplogroup
    # s_pos_list is the list of variants in the sequence
    raw_len = len(h_pos_list)
    if start:
        h_pos_list = [x for x in h_pos_list if x.start >= start]
        #print len(h_pos_list)
    if end:
        h_pos_list = [x for x in h_pos_list if x.start <= end]
        #print len(h_pos_list)
    # cerca se le mutazioni sono presenti nella lista
    matched = sum(1 for x in h_pos_list if x in s_pos_list)
    missing_haplo_pos = sorted([x for x in h_pos_list if x not in s_pos_list],
                               key=lambda x: x.start)
    # cerca poi se NON ci sono le retromutazioni nella lista

    # ------accrocchio da modificare
    # eliminare dal conto le posizioni retromutate
    # bisogna tener conto delle posizioni Retromutated
    # retrom = [y for y in h_pos_list if isinstance(y, (datatypes.Retromutation, datatypes.Retromutated))]
    # mut_pos = [x.start for x in s_pos_list]
    # matched += sum(1 for x in retrom if x.start not in mut_pos)
    # --------------

    return matched, raw_len, len(h_pos_list), missing_haplo_pos


def compare_mutations_regions(h_pos_list, s_pos_list, regions=None):
    """
    Finds mutations in multiple regions (e.g. for contigs)
    regions are encoded like [(1,10000), (11000,15000)]
    """
    # h_pos_list is the list of variants defining a haplogroup
    # s_pos_list is the list of variants in the sequence
    raw_len = len(h_pos_list)
    h_pos_list_checked = []
    if regions:
        for region in regions:
            region_start, region_end = region[0], region[1]
            h_pos_list_checked.extend([x for x in h_pos_list
                                       if region_start <= x <= region_end])
    else:
        h_pos_list_checked = h_pos_list
    matched = sum(1 for x in h_pos_list_checked if x in s_pos_list)
    missing_haplo_pos = sorted([x for x in h_pos_list_checked if x not in s_pos_list],
                               key=lambda x: x.start)
    return matched, raw_len, len(h_pos_list_checked), missing_haplo_pos


def align_sequences(rif, inc):  # muscle
    muscle = subprocess.Popen(['./muscle -quiet'], shell=True,
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              universal_newlines=True)
    alg = muscle.communicate(">%s\n%s\n>%s\n%s\n" % (rif.name, rif.seq,
                                                     inc.name, inc.seq))[0]
    # print alg
    alg_split = alg.split('>')[1:]
    rif_alg = ''.join(alg_split[0].split()[1:]).upper()
    inc_alg = ''.join(alg_split[1].split()[1:]).upper()
    return (alg, Sequence(rif.name, rif_alg), Sequence(inc.name, inc_alg))

# Alcune Sequenze non Vengono correttamente allineate da Muscle, sarebbe da implementare anche mafft?
# def align_sequences_mafft(rif, inc):
#     mafft = subprocess.Popen(['mafft'], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
#     alg = muscle.communicate(">%s\n%s\n>%s\n%s\n" % (rif.name, rif.seq, inc.name, inc.seq))[0]
#     alg_split = alg.split('>')[1:]
#     rif_alg = ''.join(alg_split[0].split()[1:]).upper()
#     inc_alg = ''.join(alg_split[1].split()[1:]).upper()
#     return alg, datatypes.Sequence(rif.name, rif_alg), datatypes.Sequence(inc.name, inc_alg)


# TODO: check if these are class or instance variables
class SequenceDiff(object):
    diff_list = None
    diff_list_raw = None
    rif = None
    obj = None
    start = None
    end = None
    gen_snps = staticmethod(get_snps)
    gen_alg = staticmethod(align_sequences)
    # eliminare Deletion(3107), Deletion(523) dall'analisi
    # ignore_position = None
    ignore_position = []

    def __init__(self, mit=True):
        if mit:
            self.ignore_position.extend([Deletion("3107d"),
                                         Deletion("523-524d")])

    def gen_diff(self, rif=None, obj=None):
        # temporaneo, per i casi in cui viene dato a parte l'allineamento
        if rif:
            alg, self.rif, self.obj = self.gen_alg(rif, obj)
            self.alg = alg
        self.diff_list = self.gen_snps(self.rif.seq, self.obj.seq)
        # print self.diff_list
        if self.ignore_position:
            for i in self.ignore_position:
                if i in self.diff_list:
                    del self.diff_list[self.diff_list.index(i)]
                    print('_'*30)
                    print("**** Deleting SNP: %s" % i.pprint())

    def find_segment(self):
        if not self.rif is None:
            b_idx = min(self.obj.seq.find(x) for x in DNA)
            # funge fino a quando il blocco non contiene dei gap
            #e_idx = self.obj.seq.find('-', b_idx)
            tmp = ''.join(reversed(self.obj.seq))
            e_idx = len(self.obj.seq) - min(tmp.find(x) for x in DNA)
            b_idx = len(self.rif.seq[:b_idx].replace('-', ''))
            e_idx = len(self.rif.seq[:e_idx].replace('-', ''))
            self.start = b_idx
            self.end = e_idx
            self.diff_list_raw = self.diff_list
            self.diff_list = [x for x in self.diff_list
                              if b_idx <= x.start <= e_idx]

    def raw_print(self, outfile=sys.stdout):
        outfile.write('%s | %s\n' % (self.rif.name, self.obj.name))
        for pos, (x, y) in enumerate(zip(self.rif.seq, self.obj.seq)):
            if x != y:
                outfile.write('% 6d | %s->%s\n' % (pos + 1, x, y))

    def pprint(self, outfile=sys.stdout):
        outfile.write('%s | %s\n' % (self.rif.name, self.obj.name))
        for pos in self.diff_list:
            outfile.write('%s\n' % pos.pprint())

    def print_alg(self, outfile=sys.stdout):
        outfile.write(self.alg)

    def __repr__(self):
        return self.pprint()


class HaplogroupStats(object):
    def __init__(self, name='', count=0, total=0, raw_total=0):
        self.name = name
        self.count = count
        self.total = total
        self.raw_total = raw_total
        if total > 0:
            self.percentage = count / float(total) * 100
        else:
            self.percentage = 0

    def __eq__(self, other):
        return True if self.percentage == other.percentage else False

    def __lt__(self, other):
        return True if self.percentage < other.percentage else False

    def __gt__(self, other):
        return not self < other

    def __str__(self):
        return "%s: %.2f%%" % (self.name, self.percentage)

    def __repr__(self):
        return str(self)


class Classify(object):
    """
    Problema coi segmenti. Le delezioni che cadono in quel punto non vengono
    considerate (non di interesse)
    """
    # TODO
    # controllare il problema dell'inizializzazione di questo tipo di variabili in una class
    # stat_list = []
    # haplo_stats = {}
    # MEMO
    # matched: Nph
    # tot: Nph_exp
    # raw_tot: Nph_tot
    # haplo_stats = {'haplo_name' : (Nph, Nph_exp, Nph_tot)}
    # haplogroup prediction % is calculated as matched/tot*100
    # haplogroup prediction sorting 
    def __init__(self):
        self.stat_list = []
        self.haplo_stats = {}
        self.missing_haplo_pos = {}

    def classify_haplogroup(self, haplo_name, haplo_list, seq_diff):
        """Classifica un SeqDiff per un dato aplogruppo"""
        matched, raw_tot, tot, missing_haplo_pos = compare_mutations(
            haplo_list, seq_diff.diff_list, seq_diff.start, seq_diff.end
        )
        if matched > 0:
            self.haplo_stats[haplo_name] = (matched, tot, raw_tot)
            self.stat_list.append(HaplogroupStats(haplo_name, matched, tot, raw_tot))
            self.missing_haplo_pos[haplo_name] = missing_haplo_pos

    def classify_by_tree(self, haplo_tree, seq_diff, filt=True):
        """Classifica un SeqDiff usando un albero creato in precedenza"""
        for haplo_name in haplo_tree:
            if filt:
                h_pos_list = haplo_tree.get_filtered_positions(haplo_name)
            else:
                h_pos_list = haplo_tree.get_branch_positions(haplo_name)
            self.classify_haplogroup(haplo_name, h_pos_list, seq_diff)
        self.seq_diff = seq_diff

    def classify90(self):
        """Generate haplo_stats for which Nph/Nph_exp > 0.89, sorted for decreasing P_Hg"""
        self.haplo_stats90 = dict((key, value) for key, value in self.haplo_stats.iteritems()
                                  if value[0]/float(value[1]) > 0.89)
        self.haplo_stats90 = OrderedDict(sorted(self.haplo_stats90.items(),
                                                key=lambda x: x[1][0]/float(x[1][1]),
                                                reverse=True))
        return self.haplo_stats90

    def best_predictions(self):
        self.haplo_stats90_list = sorted(self.haplo_stats90.items(),
                                         key=lambda x: x[1][0], reverse=True)

    def pprint(self, outfile=sys.stdout):
        """Output della classificazione, è un csv"""
        outfile.write("Sequence Name,Predicted Haplogroup,N,Nph,Nph_tot,Nph_exp,P_Hg\n")
        seq_name = self.seq_diff.obj.name
        seq_snps = len(self.seq_diff.diff_list)
        for haplo in self.haplo_stats:
            count, total, raw_tot = self.haplo_stats[haplo]
            if total > 0:
                outfile.write("%s,%s,%d,%d,%d,%d,%.3f\n" % (seq_name, haplo, seq_snps,
                                                            count, raw_tot, total,
                                                            count/float(total)*100))

    def pprint_sorted(self, outfile=sys.stdout):
        """Output della classificazione, è un csv"""
        outfile.write("Sequence Name,Predicted Haplogroup,N,Nph,Nph_tot,Nph_exp,P_Hg,Missing sites\n")
        seq_name = self.seq_diff.obj.name
        seq_snps = len(self.seq_diff.diff_list)
        for haplo in self.haplo_stats_sorted.keys():
            count, total, raw_tot = self.haplo_stats_sorted[haplo]
            missing_haplo_pos_tostring = ';'.join([i.__str__()
                                                   for i in self.missing_haplo_pos[haplo]])
            outfile.write("%s,%s,%d,%d,%d,%d,%.3f,%s\n" % (seq_name, haplo, seq_snps,
                                                           count, raw_tot, total,
                                                           count/float(total)*100,
                                                           missing_haplo_pos_tostring))

    def get_genome_state(self):
        """
        Check if the genome is complete.
        The check is performed on alignment start/end,
        but could be performed also by checking Nph_exp == Nph_tot
        """
        if self.seq_diff.start == 0 and self.seq_diff.end == 16569:
            self.genome = "complete"
        else:
            self.genome = "incomplete"
        return self.genome

    def prediction_sorting(self):
        s = self.get_genome_state()
        if s == "complete":
            t = "complete"
            # apply sorting P_Hg decreasing, then Nph decreasing
            self.haplo_stats_sorted = OrderedDict(sorted(self.haplo_stats.items(),
                                                         key=lambda x: x[1][0], reverse=True))
            self.haplo_stats_sorted = OrderedDict(sorted(self.haplo_stats.items(),
                                                         key=lambda x: x[1][0]/float(x[1][1]),
                                                         reverse=True))
            # check if there are two haplogroups with max values of Nph and P_Hg
            # TODO: this is 1 metre long, better to refactor
            self.haplo_best = OrderedDict(
                filter(
                    lambda x: (x[1][0], x[1][0]/float(x[1][1])) == (self.haplo_stats_sorted.items()[0][1][0],
                                                                    self.haplo_stats_sorted.items()[0][1][0] / float(self.haplo_stats_sorted.items()[0][1][1])),
                    self.haplo_stats_sorted.items()))
        else:
            t = "incomplete"
            # apply sorting: take res with P_Hg > 90,then the one(s) with
            # highest Nph is/are the best
            if self.classify90():
                # print "got 90"
                self.haplo_stats_sorted = self.classify90()
                self.haplo_stats_sorted = OrderedDict(sorted(self.haplo_stats_sorted.items(),
                                                             key=lambda x: x[1][0],
                                                             reverse=True))
                # TODO: see above
                self.haplo_best = OrderedDict(filter(lambda x: (x[1][0],x[1][0]/float(x[1][1])) == (self.haplo_stats_sorted.items()[0][1][0],self.haplo_stats_sorted.items()[0][1][0]/float(self.haplo_stats_sorted.items()[0][1][1])), self.haplo_stats_sorted.items()))

            # cases with no P_Hg >= 0.9, just take the haplogroup(s) with the highest P_Hg
            else:
                self.haplo_stats_sorted = OrderedDict(sorted(self.haplo_stats.items(),
                                                             key=lambda x: x[1][0]/float(x[1][1]),
                                                             reverse=True))
                self.haplo_best = OrderedDict(filter(lambda x: x[1][0]/float(x[1][1]), self.haplo_stats_sorted.items()))
        return self.haplo_stats_sorted, self.haplo_best

    def get_marges(self, marge_dict):
        print("hidden best is", self.haplo_best)
        best_list = [i for i in self.haplo_best]
        print("best_list is", best_list)
        # TODO: I guess the following is parse_mhcs.which_mhcs_lite()
        #   instead of parse_marge.which_marge_lite()
        self.marges = set([which_mhcs_lite(i, marge_dict)
                           for i in best_list])
        print("marges are ", self.marges)
        return self.marges
