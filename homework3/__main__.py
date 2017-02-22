#!/usr/bin/env python3

"""
Align sequences using a given substitution matrix

Usage:
    align.py <substitution_matrix> <sequence_pairs> [options]
    align.py gaps
    align.py thresholds
    align.py compare [options]
    align.py optimize


Arguments:
    <substitution_matrix>
        Name of the substitution matrix to use for the alignments

    <sequence_pairs>
        Name of file containing space delimited pairs of sequences to align

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


Options:
    -n --normalize
        Normalize the raw scores from the alignment by the length of the
        shorter of the two sequences

    -o --output <path>
        Save alignment output to a file named <path>

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

    record_thresholds.to_csv('Thresholds.csv')

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
        if normalize == True:
            seq_align = align(normalize=True)
            print('\n\n\n\n\n')
            print("NORMALIZING SCORES")
            print(seq_align.normalize)
        else:
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
        if normalize == True:
            seq_align = align(normalize=True)
        else:
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
    matrix_false_pos.to_csv("Compare_matrices.csv")


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

        print(score_dict[matrix]['tp'])
        print(score_dict[matrix]['fp'])
        print(ruler_set)

        ruler = sorted(list(ruler_set), reverse=True)

        # Get counts of values above threshold in tp and fp lists
        tp_threshold_count = len([element for element in score_dict[matrix]['tp'] if element > score_dict[matrix]['threshold']])
        fp_threshold_count = len([element for element in score_dict[matrix]['fp'] if element > score_dict[matrix]['threshold']])

        x = [] # False positive
        y = [] # True positive

        # Count and append hits to x-axis and y-axis lists
        for tick in ruler:
            x.append(len([element for element in score_dict[matrix]['fp'] if element >= tick and element >= score_dict[matrix]['threshold']])/fp_threshold_count)
            y.append(len([element for element in score_dict[matrix]['tp'] if element >= tick and element >= score_dict[matrix]['threshold']])/tp_threshold_count)

        plt.plot(x, y, color=colors[matrix], lw=lw, label=matrix)

    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic - Raw')
    plt.legend(loc="lower right")
    plt.show()
    fig.savefig('Matrix_ROC-Raw.pdf', dpi=300)


def run_alignments(seq_align, save_output=None):
    """
    Core code to run alignments between two sequences

    Parameters
    ----------
    seq_align

    """
    if save_output:
        output = open(save_output, 'w+')
        write_counter = 1

    for working_pair in seq_align.working_pairs:
        seq_align.seq_A = SeqIO.read(working_pair.split()[0], 'fasta').upper()
        seq_align.seq_B = SeqIO.read(working_pair.split()[1], 'fasta').upper()

        print(seq_align.seq_A.seq)
        print(seq_align.seq_B.seq)

        # Do work
        seq_align.initialize_scoring_matrix()
        seq_align.fill_in_matrix(working_pair)
        seq_A_string, bar_things, seq_B_string = seq_align.traceback_II()

        if save_output:
            output.write('>{}\n'.format(write_counter))
            write_counter += 1
            output.write(seq_A_string + '\n')
            output.write(bar_things + '\n')
            output.write(seq_B_string + '\n')


def _crappy_parser(file_name):

        align = open(file_name)
        seq_list = []
        temp_list = []

        for line in align:
            if line[0] != '>':
                temp_list.append(line.strip())

            if len(temp_list) == 3:
                seq_list.append(temp_list)
                temp_list = []

        print("\n{} imported".format(file_name))

        return seq_list


def matrix_optimization():
    """
    Optimize the matrix... somehow...
    Objective function: sum of TPR for FPRs 0, 0.1, 0.2, 0.3

    SIMULATED ANNEALING MONTE CARLO WITH METROPOLIS HASTINGS CRITERIA
    Sample weighted for AA frequency in alignment sequences
    Select random AA pairs, change value in matrix
    Calculate objective function value
        If current < previous, accept new matrix
        Else, calculate boltzmann probability (e^-(d(score / kT))), accept if > ( 0 > randint > 1 )
    Increment temperature every 100 or so iterations

    This is going to be super slow...
    """

    # Set substitution matrix to PAM100
    working_substitution_matrix = pd.read_table(open('./BLOSUM50'), delim_whitespace=True, header=6)
    possible_AA = working_substitution_matrix.columns.values
    working_substitution_matrix = working_substitution_matrix.set_index(possible_AA)

    # Set gap penalties
    gap_opening = -6
    gap_extension = -4

    # Import saved positive and negative alignments and save to list of lists
    # [Seq_A, bar_things, Seq_B]
    neg_seq_list = _crappy_parser('Alignments-neg.txt')
    pos_seq_list = _crappy_parser('Alignments-pos.txt')

    # Calculate AA frequencies in Positive alignments
    all_pos_align_sequences = [element[0] + element[2] for element in pos_seq_list]
    one_big_pos_sequence = re.sub('[^ARNDCQEGHILKMFPSTWYVBZX*]', '', ''.join(all_pos_align_sequences))

    mat_dict = {'A': 0,
                'R': 1,
                'N': 2,
                'D': 3,
                'C': 4,
                'Q': 5,
                'E': 6,
                'G': 7,
                'H': 8,
                'I': 9,
                'L': 10,
                'K': 11,
                'M': 12,
                'F': 13,
                'P': 14,
                'S': 15,
                'T': 16,
                'W': 17,
                'Y': 18,
                'V': 19,
                'B': 20,
                'Z': 21,
                'X': 22,
                '*': 23
                }

    # MC!
    accepted_matrix = None
    current_matrix = working_substitution_matrix.values
    accepted_TPR = 0
    temperatures = [1, 0.8, 0.6, 0.4, 0.2, 0.1]

    import random

    # Run MC
    for temperature in temperatures:
        print("Current Temp: {}".format(temperature))
        for increment in range(100):
            # Generate new matrix
            first_AA = mat_dict[random.choice(one_big_pos_sequence)]
            second_AA = mat_dict[random.choice(one_big_pos_sequence)]

            # Maintain symmetry
            current_matrix[first_AA, second_AA] = current_matrix[first_AA, second_AA] + random.uniform(-1, 1) * temperature
            current_matrix[second_AA, first_AA] = current_matrix[first_AA, second_AA] + random.uniform(-1, 1) * temperature

            # Get FP thresholds from saved alignment
            FP_scores = []
            current_score = []
            extending = False

            for alignment in neg_seq_list:

                # print(alignment[0])
                # print(alignment[1])
                # print(alignment[2])

                for seq_A, seq_B in zip(alignment[0], alignment[2]):
                    if seq_A == '-' or seq_B == '-':
                        if extending == False:
                            current_score.append(gap_opening)
                            extending = True
                        else:
                            current_score.append(gap_extension)
                            pass
                    else:
                        current_score.append(current_matrix[mat_dict[seq_A], mat_dict[seq_B]])
                        extending = False

                print(current_score)
                print(sum(current_score))

                FP_scores.append(sum(current_score))
                current_score = []

            # Calculate TP scores and get TPR for each FP threshold
            FP_thresholds = [sorted(FP_scores)[49], sorted(FP_scores)[44], sorted(FP_scores)[39], sorted(FP_scores)[34]]
            print('FP_Thresholds: {}'.format(FP_thresholds))

            TP_scores = []
            current_score = []
            extending = False

            for alignment in pos_seq_list:
                for seq_A, seq_B in zip(alignment[0], alignment[2]):
                    if seq_A == '-' or seq_B == '-':
                        if extending == False:
                            current_score.append(gap_opening)
                            extending = True
                        else:
                            current_score.append(gap_extension)
                            pass
                    else:
                        current_score.append(current_matrix[mat_dict[seq_A], mat_dict[seq_B]])
                        extending = False

                print(current_score)
                print(sum(current_score), '\n')

                TP_scores.append(sum(current_score))
                current_score = []

            sum_TPR = sum([len([element for element in TP_scores if element > threshold])/50 for threshold in FP_thresholds])

            print([len([element for element in TP_scores if element > threshold])/50 for threshold in FP_thresholds])
            print(TP_scores)
            print('sum_TPR: {}'.format(sum_TPR))

            # Reject/Accept
            if sum_TPR > accepted_TPR:
                accepted_TPR = sum_TPR
                accepted_matrix = current_matrix
            else:
                pass_limit = random.uniform(0, 1)
                boltz = np.exp(-((sum_TPR - accepted_TPR)/temperature))
                if boltz > pass_limit:
                    accepted_TPR = sum_TPR
                    accepted_matrix = current_matrix

    optimized_matrix = np.savetxt("Optimized_matrix.csv", accepted_matrix, delimiter=",")


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
        # Initialize variables and stuff

        substitution_matrices = {'BLOSUM50': 6,
                                 'BLOSUM62': 6,
                                 'MATIO': 2,
                                 'PAM100': 9,
                                 'PAM250': 9
                                 }

        seq_align.substitution_matrix = pd.read_table(open('./{}'.format(args['<substitution_matrix>'])),
                                                      delim_whitespace=True,
                                                      header=substitution_matrices[args['<substitution_matrix>']]
                                                      )

        seq_align.substitution_matrix = seq_align.substitution_matrix.set_index(seq_align.substitution_matrix.columns.values)

        seq_align.substitution_matrix.to_csv('{}.csv'.format(args['<substitution_matrix>']))

    # Import sequences to align
    if args['<sequence_pairs>']:
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

        compare_matrices(normalize)

    if args['optimize'] == True:
        matrix_optimization()