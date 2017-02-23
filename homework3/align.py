import os
import sys
import re
import numpy as np
import pandas as pd
from Bio import SeqIO

class align(object):
    def __init__(self, normalize=False):
        self.substitution_matrix = None
        self.scoring_matrix = None
        self.traceback_matrix = None
        self.seq_A = None
        self.seq_B = None

        # Gap opening and extension related
        self.gap_opening_penalty = -7
        self.gap_extension_penalty = -3
        # Keep track of gap openings across rows (iterate over rows to fill matrix so this is easy)
        self.extending_gap = False

        self.normalize = normalize
        self.max_alignment_score = []

    def initialize_scoring_matrix(self):
        """
        Initialize a [len(seq_A + 1) x len(seq_B + 1)] scoring matrix where columns and row indicies
        are 1-indexed (0-index reserved for 0 scores). Also initialize a traceback matrix to
        keep track of decisions (0: diag, 1: left, 2: up)

        Parameters
        ----------
        seq_A - represented across COLUMNS
        seq_B - represented down ROWS

        """

        # self.scoring_matrix = pd.DataFrame(columns= np.arange(len(self.seq_A) + 1), index= np.arange(len(self.seq_B) + 1))
        self.scoring_matrix = np.full((len(self.seq_B) + 1, len(self.seq_A) + 1), -1) # [ rows x columns ]
        self.traceback_matrix = np.full((len(self.seq_B) + 1, len(self.seq_A) + 1), -1) # [ rows x columns ]


    def fill_in_matrix(self, working_pair):
        """
        Fill in the scoring matrix based on the user selected substitution matrix

        """

        print('\n', working_pair)

        # Initialize 0-indexed row and column
        for column in np.arange(len(self.seq_A) + 1):
            self.scoring_matrix[0, column] = 0
        for row in np.arange(len(self.seq_B) + 1):
            self.scoring_matrix[row, 0] = 0

        # Make a list to keep track of column gap openings
        seq_A_gap_openings = [False] * (len(self.seq_A) + 1)

        # Iterate over rows and fill in the matrix based on substitution matrix
        for row_index, row in enumerate(self.scoring_matrix):

            # Reset row gap opening bool at start of each new row
            extending_gap = False

            for column_index in np.arange(len(row)):
                if row_index > 0 and column_index > 0:

                    # Determine similarity score
                    similarity_score = self.substitution_matrix.ix[str(self.seq_B.seq[row_index - 1]), str(self.seq_A.seq[column_index - 1])]

                    # Decide whether to use gap_extension or gap_opening penalty score across rows
                    if extending_gap == False:
                        gap_penalty_left = self.gap_opening_penalty
                    else:
                        gap_penalty_left = self.gap_extension_penalty

                    # Decide whether to use gap extension or gap opening penalty score down columns
                    if seq_A_gap_openings[column_index] == False:
                        gap_penalty_up = self.gap_opening_penalty
                    else:
                        gap_penalty_up = self.gap_extension_penalty

                    # Calculate possible scores for current cell in scoring matrix
                    score_up = self.scoring_matrix[row_index - 1, column_index] + gap_penalty_up
                    score_left = self.scoring_matrix[row_index, column_index - 1] + gap_penalty_left
                    score_diag = self.scoring_matrix[row_index - 1, column_index - 1] + similarity_score

                    # Take highest final score and assign to current in scoring matrix
                    final_score = max(score_up, score_left, score_diag)

                    # Toggle gap_opening or gap_extension across rows and down columns
                    # Priority diag > left > up
                    if max(score_up, score_left, score_diag) == score_diag:
                        extending_gap = False
                        seq_A_gap_openings[column_index] = False
                        self.traceback_matrix[row_index, column_index] = 0
                    elif max(score_up, score_left, score_diag) == score_left:
                        extending_gap = True
                        seq_A_gap_openings[column_index] = False
                        self.traceback_matrix[row_index, column_index] = 1
                    else: # max(score_up, score_left, score_diag) == score_up
                        seq_A_gap_openings[column_index] = True
                        extending_gap = False
                        self.traceback_matrix[row_index, column_index] = 2

                    self.scoring_matrix[row_index, column_index] = final_score if final_score > 0 else 0

        assert len(np.argwhere(self.scoring_matrix == -1)) == 0, "The scoring matrix was not completely filled in!"

        np.savetxt("Traceback_matrix.csv", self.traceback_matrix, delimiter=",")


        names = re.split('/|\.', working_pair.strip())
        file_name = '_'.join([names[1], names[3]])
        np.savetxt("{}.csv".format(file_name), self.scoring_matrix, delimiter=",")


    def traceback_II(self):
        # Find maximum score in matrix

        row_max_indicies = np.amax(self.scoring_matrix, axis=0)
        seq_B_high_score = np.argmax(row_max_indicies, axis=0)
        seq_A_high_score = np.argmax(self.scoring_matrix, axis=0)

        column_index_high = seq_B_high_score
        row_index_high = seq_A_high_score[column_index_high]

        print("Starting Column index: {}".format(column_index_high))
        print("Starting Row index: {}".format(row_index_high))

        # Initialize variables needed for traceback
        current_score = self.scoring_matrix[row_index_high, column_index_high]

        # Record maximum scores for each alignment
        if self.normalize == True:
            shorter_length = min(len(self.seq_A.seq), len(self.seq_B.seq))
            self.max_alignment_score.append(current_score / shorter_length)

            print("SHORTER LENGTH: {}".format(shorter_length))
            print("NORMALIZED SCORE: {}".format(current_score / shorter_length))

        else:

            self.max_alignment_score.append(current_score)

        print("Current Score: {}".format(current_score))

        seq_A_alignment = []
        seq_B_alignment = []

        current_row_index = row_index_high
        current_column_index = column_index_high

        # Begin traceback, stop when (above | left | diagonal) == 0
        while self.scoring_matrix[current_row_index, current_column_index] > 0:
            direction = self.traceback_matrix[current_row_index, current_column_index]

            if direction == 0: # Diag
                seq_A_alignment.append(str(self.seq_A.seq[current_column_index - 1]))
                seq_B_alignment.append(str(self.seq_B.seq[current_row_index - 1]))
                current_row_index -= 1
                current_column_index -= 1

            elif direction == 1: # Left
                seq_A_alignment.append(str(self.seq_A.seq[current_column_index - 1]))
                seq_B_alignment.append('-')
                current_column_index -= 1

            elif direction == 2: # Up
                seq_A_alignment.append('-')
                seq_B_alignment.append(str(self.seq_B.seq[current_row_index - 1]))
                current_row_index -= 1

            else:
                print('How did you get here')

        # construct_bar_things
        bar_thing_list = []
        for seq_A_seq, seq_B_seq in zip(seq_A_alignment, seq_B_alignment):
            if seq_A_seq == seq_B_seq:
                bar_thing_list.append('|')
            elif seq_A_seq == '-' or seq_B_seq == '-':
                bar_thing_list.append(' ')
            else:
                bar_thing_list.append(':')

        print(''.join(reversed(seq_A_alignment)))
        print(''.join(reversed(bar_thing_list)))
        print(''.join(reversed(seq_B_alignment)))

        seq_A_string = ''.join(reversed(seq_A_alignment))
        bar_things = ''.join(reversed(bar_thing_list))
        seq_B_string = ''.join(reversed(seq_B_alignment))

        return seq_A_string, bar_things, seq_B_string