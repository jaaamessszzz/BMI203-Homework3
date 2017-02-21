#!/usr/bin/env python3

"""
Align sequences using a given substitution matrix

Usage:
    align.py <substitution_matrix> <sequence_pairs> [options]
    align.py gaps
    align.py thresholds
    align.py compare [options]


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

    compare
        Compare substitution matrices in terms of false positive rate. Also
        generate ROC curves


Options:
    -n --normalize
        Normalize the raw scores from the alignment by the length of the
        shorter of the two sequences

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

            # Set substitution matrix to BLOSUM50
            seq_align.substitution_matrix = pd.read_table(open('./BLOSUM50'), delim_whitespace=True, header=6)
            seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

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

        # Set substitution matrix to BLOSUM50
        seq_align.substitution_matrix = pd.read_table(open('./BLOSUM50'), delim_whitespace=True, header=6)
        seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

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

def compare_matrices(normalize):
    """
    Compare the different substitution matrices in terms of false positives

    Parameters
    ----------
    normalize

    """
    # Initialize variables and stuff
    substitution_matrices = {'BLOSUM50': 6,
                             'BLOSUM62': 6,
                             'MATIO': 2,
                             'PAM100': 9,
                             'PAM250': 9
                             }

    matrix_false_pos = pd.DataFrame(columns=['Matrix', 'False_Positive_Rate'])

    score_dict = {}

    # Loop through matricies and save score lists for true_pos and false_pos
    # Generate ROC and plot
    for sub_matrix in substitution_matrices:
        score_dict[sub_matrix] = {}

        # Find 0.7 threshold for given Matrix
        seq_align = align()
        seq_align.working_pairs = open('./Pospairs.txt')

        print('\n\n\n\n\n')
        print("OBTAINING THRESHOLD SCORE FOR {}".format(sub_matrix))
        print('\n\n\n\n\n')

        # Set substitution matrix
        seq_align.substitution_matrix = pd.read_table(open('./{}'.format(sub_matrix)), delim_whitespace=True, header=substitution_matrices[sub_matrix])
        seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

        run_alignments(seq_align)

        # Save pospair scores
        # Sort max scores for all pospairs and take 15th lowest score as 0.7 threshold
        score_dict[sub_matrix]['tp'] = sorted(seq_align.max_alignment_score)
        score_dict[sub_matrix]['threshold'] = score_dict[sub_matrix]['tp'][14]

        print('\n\n\n\n\n')
        print("{} THRESHOLD: {}".format(sub_matrix, score_dict[sub_matrix]['threshold']))

        # Find false positive rate using previously found threshold value
        seq_align = align()
        seq_align.working_pairs = open('./Negpairs.txt')

        print('\n\n\n\n\n')
        print("DETERMINING FALSE POSITIVE RATE FOR {}".format(sub_matrix))
        print('\n\n\n\n\n')

        # Set substitution matrix
        seq_align.substitution_matrix = pd.read_table(open('./{}'.format(sub_matrix)), delim_whitespace=True, header=substitution_matrices[sub_matrix])
        seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

        run_alignments(seq_align)

        # Save negpair scores
        score_dict[sub_matrix]['fp'] = sorted(seq_align.max_alignment_score)

        # Get counts for elements that get scores above threshold
        above_threshold = [element for element in seq_align.max_alignment_score if element > score_dict[sub_matrix]['threshold']]
        false_positive_rate = len(above_threshold)/50
        print("False Positive Rate: {}".format(false_positive_rate))
        new_row = pd.Series({'Matrix': sub_matrix,
                             'False_Positive_Rate': false_positive_rate
                             })

        matrix_false_pos = matrix_false_pos.append(new_row, ignore_index=True)

    generate_ROC(score_dict)
    # matrix_false_pos.to_csv("Compare_matrices.csv")


def generate_ROC(score_dict):
    """
    Plot false_pos vs. true_pos in ROC
    For a given substitution matrix:
      * Combine both true_pos and false_pos score lists in to set, list,
        reverse sort (highest to lowest)
      * Look at true_pos and false_pos lists individually and determine
        the total number of elements above score threshold. Store values in
        total_above_threshold_pos/neg
      * Iterate through sorted list of all scores (ruler)
      * For each value in ruler, look at proportion of scores below value
        considered hits
      * Append counts/total_above_threshold_pos/neg to x-axis and y-axis lists

    Plotting code lazily ripped from sklearn ROC example...
    """

    fig = plt.figure(figsize=(8, 8), facecolor='white')
    lw = 2

    colors = {'BLOSUM50': sns.xkcd_rgb["pale red"],
              'BLOSUM62': sns.xkcd_rgb["grass"],
              'MATIO': sns.xkcd_rgb["cerulean"],
              'PAM100': sns.xkcd_rgb["purplish"],
              'PAM250': sns.xkcd_rgb["golden yellow"]
              }

    for matrix in score_dict:
        # Combine true_pos and false_pos score lists
        ruler_set = set()
        for asdf in score_dict[matrix]['tp']:
            ruler_set.add(asdf)
        for asdf in score_dict[matrix]['fp']:
            ruler_set.add(asdf)


        ruler = sorted(list(ruler_set), reverse=True)
        max_value = max(ruler)

        # Get counts of values above threshold in tp and fp lists
        tp_threshold_count = len([element for element in score_dict[matrix]['tp'] if element > score_dict[matrix]['threshold']])
        fp_threshold_count = len([element for element in score_dict[matrix]['fp'] if element > score_dict[matrix]['threshold']])

        x = [] # False positive
        y = [] # True positive

        # Count and append hits to x-axis and y-axis lists
        for tick in np.arange(0, max_value + 1, 0.1): # Necessary for step function...
            x.append(len([element for element in score_dict[matrix]['fp'] if element >= tick and element >= score_dict[matrix]['threshold']])/fp_threshold_count)
            y.append(len([element for element in score_dict[matrix]['tp'] if element >= tick and element >= score_dict[matrix]['threshold']])/tp_threshold_count)

        plt.plot(x, y, color=colors[matrix], lw=lw, label=matrix)

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.show()
    fig.savefig('Matrix_ROC.pdf', dpi=300)


def run_alignments(seq_align):
    """
    Core code to run alignments between two sequences

    Parameters
    ----------
    seq_align

    """

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
    import seaborn as sns
    import matplotlib.pyplot as plt

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

    if args['compare'] == True:
        if args['--normalize']:
            normalize = True
        else:
            normalize = False

        compare_matrices(normalize)