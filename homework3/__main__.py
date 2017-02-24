#!/usr/bin/env python3

"""
Align sequences using a given substitution matrix

Usage:
    homework3 align <substitution_matrix> <sequence_pairs> [options]
    homework3 gaps
    homework3 thresholds
    homework3 compare [options]
    homework3 optimize <matrix_to_optimize>


Arguments:
    <substitution_matrix>
        Name of the substitution matrix to use for the alignments

    <sequence_pairs>
        Name of file containing space delimited pairs of sequences to align

    align
        Run alignment with specified substitution matrix and sequence pairs

    thresholds
        Run routine to determine the 0.7 score threshold for each gap opening/
        extension penalty combination

    gaps
        Run routine to determine optimal gap penalties for the BLOSUM50 matrix
        using the previously determined threshold scores

    compare
        Compare substitution matrices in terms of false positive rate. Also
        generate ROC curves

    optimize
        Run algorithm to optimize scoring matrix...

    <matrix_to_optimize>
        Name of the matrix to run optimization on


Options:
    -n --normalize
        Normalize the raw scores from the alignment by the length of the
        shorter of the two sequences

    -o --output <path>
        Save alignment output to a file named <path>

    -c --compare_optimized <matrix>
        Compare 1) default matrix, 2) optimized scoring matrix against default
        matrix alignments, and 3) optimized scoring matrix against optimized
        alignments

"""

if __name__ == '__main__':
    from .align import align
    from .util import determine_thresholds, determine_gap_penalties, run_alignments, compare_matrices, compare_optimized, matrix_optimization

    import docopt
    import re
    import collections
    import os
    import sys
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
    import seaborn as sns
    import matplotlib.pyplot as plt

    args = docopt.docopt(__doc__)
    seq_align = align()

    # Set substitution matrix
    if args['align']:
        # Initialize variables and stuff

        substitution_matrices = {'BLOSUM50': 6,
                                 'BLOSUM62': 6,
                                 'BLOSUM62-Optimized': 0,
                                 'MATIO': 2,
                                 'MATIO-Optimized': 0,
                                 'PAM100': 9,
                                 'PAM100-Optimized': 0,
                                 'PAM250': 9
                                 }

        seq_align.substitution_matrix = pd.read_table(open('./{}'.format(args['<substitution_matrix>'])),
                                                      delim_whitespace=True,
                                                      header=substitution_matrices[args['<substitution_matrix>']]
                                                      )

        seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

        seq_align.substitution_matrix.to_csv('{}.csv'.format(args['<substitution_matrix>']))

        seq_align.working_pairs = open(args['<sequence_pairs>'])
        run_alignments(seq_align, args['--output'])

    if args['thresholds'] == True:
        determine_thresholds()

    if args['gaps'] == True:
        determine_gap_penalties()

    if args['compare'] == True:
        if args['--normalize']:
            normalize = True
        else:
            normalize = False

        if args['--compare_optimized']:
            compare_optimized(args['--compare_optimized'])
        else:
            compare_matrices(normalize)

    if args['optimize'] == True:
        matrix_optimization(args['<matrix_to_optimize>'])