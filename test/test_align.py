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

    for test_pair in test_pairs:
        print(test_pair)
        seq_align.working_pairs = [test_pair]
        A, bar, B = run_alignments(seq_align)

        if test_pair == "test/test_1.fa  test/test_2.fa":
            assert A == "RRRRRAAAAA", (A, type(A))
            assert bar == "||||||||||", (bar, type(bar))
            assert B == "RRRRRAAAAA", (B, type(B))

        elif test_pair == "test/test_3.fa  test/test_4.fa":
            assert A == None, (A, type(A))
            assert bar == None, (bar, type(bar))
            assert B == None, (B, type(B))

        elif test_pair == "test/test_1.fa  test/test_3.fa":
            assert A == "WWWWWWWWWWWWWWW", (A, type(A))
            assert bar == "|||||||||||||||", (bar, type(A))
            assert B == "WWWWWWWWWWWWWWW", (B, type(A))

        elif test_pair == "test/test_2.fa  test/test_4.fa":
            assert A == "AAAAAAAAAAAAAAA", (A, type(A))
            assert bar == "|||||||||||||||", (bar, type(A))
            assert B == "AAAAAAAAAAAAAAA", (B, type(A))

        elif test_pair == "test/test_5.fa  test/test_6.fa":
            assert A == "ASDFASDFASDF--DFASD--SDFASDFASDFASDF", (A, type(A))
            assert bar == "||:|||:|||:|  :|||:  |:|||:|||:|||:|", (bar, type(A))
            assert B == "ASHFASHFASHFASHFASHFASHFASHFASHFASHF", (B, type(A))
            
        else:
            assert False==True, "How did you get here"

