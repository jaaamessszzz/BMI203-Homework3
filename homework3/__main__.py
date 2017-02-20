#!/usr/bin/env python3

"""
Align sequences using a given substitution matrix

Usage:
    align.py <substitution_matrix> <sequence_pairs> [options]
    align.py gaps
    align.py thresholds

Arguments:
    <substitution_matrix>
        Name of the substitution matrix to use for the alignments

    <sequence_pairs>
        Name of file containing space delimited pairs of sequences to align

    thresholds
        Run routine to determine the 0.7 score threshold for each gap opening/
        extension penalt combination

    gaps
        Run routine to determine optimal gap penalties for the BLOSUM50 matrix
        using the previously determined threshold scores


Options:
    -o --output <path>
        Output the alignments to a file named <path>
        NOT YET IMPLEMENTED

"""

def determine_thresholds():
    """
    Run routine to determine the positive pair thresholds for all given opening and extension
    penalty scores with the BLOSUM50 matrix. The 0.7 threshold score for each opening/extension
    score is saved to a pandas dataframe and exported to a csv so I only need to run this once...
    """
    # Record in... a dataframe? Everything in dataframes.
    record_thresholds = pd.DataFrame(columns=['Gap_opening', 'Gap_extension', 'Threshold_Score'])

    for gap_opening in range(1, 21, 1):
        for gap_extension in range(1,6,1):
            # I was having a problem where the run_alignments() function was only being run once in the for loop...
            # Not entirely sure why that was, but reinitializing seq_align every time is a solution
            # Not a good solution, but a solution for the time being...

            seq_align = align()

            # Import pospairs and try all gap opening and extention penalties
            seq_align.working_pairs = open('./Pospairs.txt')

            # Set gap opening and extension penalties
            seq_align.gap_opening_penalty = -(gap_opening)
            seq_align.gap_extension_penalty = -(gap_extension)

            print('\n\n\n\n\n')
            print("TESTING THE FOLLOWING GAP OPENING AND EXTENSION PENALTIES:")
            print("GAP OPENING: {}".format(seq_align.gap_opening_penalty))
            print("GAP EXTENSION: {}".format(seq_align.gap_extension_penalty))
            print('\n\n\n\n\n')

            run_alignments(seq_align)

            # Sort max scores for all pospairs and take 15th lowest score as 0.7 threshold
            threshold_score = sorted(seq_align.max_alignment_score)[14]
            print(sorted(seq_align.max_alignment_score))
            new_row = pd.Series({'Gap_opening': gap_opening,
                                 'Gap_extension': gap_extension,
                                 'Threshold_Score': threshold_score
                                 })

            record_thresholds = record_thresholds.append(new_row, ignore_index=True)

    record_thresholds.to_csv('Thresholds-asdf.csv')

def determine_gap_penalties():
    """
    Use the threshold values in Thresholds.csv to find the number of false positives in Negpairs.txt
    using the associated gap penalty values.
    Returns
    -------

    """
    gap_thresholds = pd.read_csv('Thresholds.csv')
    for index, row in gap_thresholds.iterrows():
        seq_align = align()

        # Import pospairs and try all gap opening and extention penalties
        seq_align.working_pairs = open('./Negpairs.txt')

        # Set gap opening and extension penalties
        seq_align.gap_opening_penalty = -(row['Gap_opening'])
        seq_align.gap_extension_penalty = -(row['Gap_extension'])

        print('\n\n\n\n\n')
        print("TESTING THE FOLLOWING GAP OPENING AND EXTENSION PENALTIES:")
        print("GAP OPENING: {}".format(seq_align.gap_opening_penalty))
        print("GAP EXTENSION: {}".format(seq_align.gap_extension_penalty))
        print('\n\n\n\n\n')

        run_alignments(seq_align)

        # Get counts for elements that get scores above threshold
        above_threshold = [element for element in seq_align.max_alignment_score if element > row['Threshold_Score']]

        print("CURRENT THRESHOLD: {}".format(row['Threshold_Score']))
        print(seq_align.max_alignment_score)
        print(above_threshold)

        false_positive_rate = len(above_threshold)/50
        gap_thresholds.ix[index, 'False_Positive_Rate'] = false_positive_rate

    gap_thresholds.to_csv('False_Positives.csv')

def run_alignments(seq_align):
    # Set substitution matrix to BLOSUM50
    seq_align.substitution_matrix = pd.read_table(open('./BLOSUM50'), delim_whitespace=True, header=6)
    seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

    for working_pair in seq_align.working_pairs:
        seq_align.seq_A = SeqIO.read(working_pair.split()[0], 'fasta').upper()
        seq_align.seq_B = SeqIO.read(working_pair.split()[1], 'fasta').upper()

        print(seq_align.seq_A.seq)
        print(seq_align.seq_B.seq)

        # Do work
        seq_align.initialize_scoring_matrix()
        seq_align.fill_in_matrix()
        seq_align.traceback()


if __name__ == '__main__':
    from .align import align
    import docopt
    import re
    import collections
    import os
    import sys
    import numpy as np
    import pandas as pd
    from Bio import SeqIO

    args = docopt.docopt(__doc__)
    seq_align = align()

    # Set substitution matrix
    if args['<substitution_matrix>']:
        seq_align.substitution_matrix = pd.read_table(open('./{}'.format(args['<substitution_matrix>'])), delim_whitespace=True, header=6)
        seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

        seq_align.substitution_matrix.to_csv('{}.csv'.format(args['<substitution_matrix>']))

    # Import sequences to align
    if args['<sequence_pairs>']:
        seq_align.working_pairs = open(args['<sequence_pairs>'])
        run_alignments(seq_align)

    if args['thresholds'] == True:
        determine_thresholds()

    if args['gaps'] == True:
        determine_gap_penalties()