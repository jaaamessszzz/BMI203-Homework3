"""
I'm only testing the alignment function... all the other ones are all to specific to test easily...
"""

from homework3.align import align
from homework3.util import determine_thresholds, determine_gap_penalties, run_alignments, compare_matrices, compare_optimized, matrix_optimization
import pandas as pd

def test_align():
    seq_align = align()

    # Import pospairs and try all gap opening and extention penalties
    test_pairs = open('./Testpairs.txt')

    # Set substitution matrix to BLOSUM50
    seq_align.substitution_matrix = pd.read_table(open('./BLOSUM50'), delim_whitespace=True, header=6)
    seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

    for index, test_pair in enumerate(test_pairs):
        seq_align.working_pairs = [test_pair]
        A, bar, B = run_alignments(seq_align)

        if index == 0:
            assert A == "RRRRRAAAAA", (A, type(A))
            assert bar == "||||||||||", (bar, type(bar))
            assert B == "RRRRRAAAAA", (B, type(B))

        elif index == 1:
            assert A == '', (A, type(A))
            assert bar == '', (bar, type(bar))
            assert B == '', (B, type(B))


        elif index == 2:
            assert A == "WWWWWWWWWWWWWWW", (A, type(A))
            assert bar == "|||||||||||||||", (bar, type(bar))
            assert B == "WWWWWWWWWWWWWWW", (B, type(B))

        elif index == 3:
            assert A == "AAAAAAAAAAAAAAA", (A, type(A))
            assert bar == "|||||||||||||||", (bar, type(bar))
            assert B == "AAAAAAAAAAAAAAA", (B, type(B))

        elif index == 4:
            assert A == "ASDFASDFASDF--DFASD--SDFASDFASDFASDF", (A, type(A))
            assert bar == "||:|||:|||:|  :|||:  |:|||:|||:|||:|", (bar, type(bar))
            assert B == "ASHFASHFASHFASHFASHFASHFASHFASHFASHF", (B, type(B))

        else:
            assert False == True, (test_pair, type(test_pair))

