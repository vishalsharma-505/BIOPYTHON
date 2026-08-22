from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
aligner = PairwiseAligner()

seq1 = Seq("TATAAT")
seq2 = Seq("ATGAATA")
aligner.mode = "local"
aligner.match_score = 5
aligner.open_gap_score = -1
aligner.extend_gap_score = -2




alignments = aligner.align(seq1,seq2)

for alignment in alignments :
    print(alignment)
