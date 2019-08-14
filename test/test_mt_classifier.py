from modules.mt_classifier import main_mt_hpred
import shutil, os
from modules.classifier import consts, NGclassify, datatypes
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def test_gen_diff():
    # obj = SeqIO.index("test/data/test.fasta", 'fasta')
    # ID = obj[list(obj.keys())[0]].id
    rif = SeqRecord(Seq(consts.RCRS), id = 'RSRS', name = 'RSRS')
    seqtest_handle = SeqIO.index("test/data/test.fasta", 'fasta')
    seqtest = seqtest_handle[list(seqtest_handle.keys())[0]]
    print(type(seqtest))
    seq_diff = NGclassify.SequenceDiff()
    seq_diff.gen_diff(muscle_exe = shutil.which('muscle'), rif = rif, obj = seqtest)
    #assert(seq_diff.diff_list == [datatypes.Transversion("36C"), datatypes.Transversion("321A"), datatypes.Transversion("4456G"), datatypes.Transversion("8243C")])
    assert(seq_diff.diff_list == [datatypes.Transition(146), datatypes.Transition(195), datatypes.Transition(235), datatypes.Transition(247), datatypes.Insertion("310.C"), datatypes.Transition(471), datatypes.Transition(663), datatypes.Transition(769), datatypes.Transversion("825T"), datatypes.Transition(1018), datatypes.Transition(1442), datatypes.Transition(1736), datatypes.Transition(2758), datatypes.Transition(2885), datatypes.Transition(3594), datatypes.Transition(4104), datatypes.Transition(4248), datatypes.Transition(4312), datatypes.Transition(4824), datatypes.Transition(6266), datatypes.Transition(7146), datatypes.Transition(7256), datatypes.Transition(7521), datatypes.Transition(8468), datatypes.Transition(8655), datatypes.Transition(8701), datatypes.Transition(8794), datatypes.Transition(9540), datatypes.Transition(10398), datatypes.Transition(10664), datatypes.Transition(10688), datatypes.Transition(10810), datatypes.Transition(10873), datatypes.Transition(10915), datatypes.Transition(11914), datatypes.Transition(13105), datatypes.Transition(13276), datatypes.Transition(13506), datatypes.Transition(13650), datatypes.Transition(16129), datatypes.Transition(16187), datatypes.Transition(16189), datatypes.Transition(16230), datatypes.Transition(16278), datatypes.Transition(16290), datatypes.Transition(16311), datatypes.Transition(16319), datatypes.Transition(16362), datatypes.Transition(16519)])

# def test_main_mt_hpred():
#     contig_file = "test/data/test.fasta"
#     muscle_exe = shutil.which('muscle')
#     basename = "test_out"
#     best_results_file = 'test_best_results.csv'
#     data_file = "data/classifier"
#     if not os.path.exists("test/data/output"):
#         os.makedirs("test/data/output")
#     main_mt_hpred(contig_file = contig_file, muscle_exe = muscle_exe, basename = basename, best_results_file = best_results_file, data_file = data_file)

