import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

class align(object):
    def __init__(self, normalize=False):
        self.substitution_matrix = None
        self.scoring_matrix = pd.DataFrame()
        self.seq_A = None
        self.seq_B = None

        self.gap_opening_penalty = -9
        self.gap_extension_penalty = -4
        self.extending_gap = False

        self.normalize = normalize
        self.max_alignment_score = []

    def initialize_scoring_matrix(self):
        """
        Initialize a [len(seq_A) x len(seq_B)] scoring matrix where columns and row indicies
        are 1-indexed (0-index reserved for 0 scores).

        Parameters
        ----------
        seq_A - represented across COLUMNS
        seq_B - represented down ROWS

        """

        # self.scoring_matrix = pd.DataFrame(columns= np.arange(len(self.seq_A) + 1), index= np.arange(len(self.seq_B) + 1))
        self.scoring_matrix = np.full((len(self.seq_B) + 1, len(self.seq_A) + 1), -1) # [ rows x columns ]

    def fill_in_matrix(self):
        """
        Fill in the scoring matrix based on the user selected substitution matrix

        """

        # Initialize 0-indexed row and column
        for column in np.arange(len(self.seq_A) + 1):
            self.scoring_matrix[0, column] = 0
        for row in np.arange(len(self.seq_B) + 1):
            self.scoring_matrix[row, 0] = 0

        # Iterate over rows and fill in the matrix based on substitution matrix
        for row_index, row in enumerate(self.scoring_matrix):
            for column_index in np.arange(len(row)):
                if row_index > 0 and column_index > 0:

                    # Determine similarity score
                    similarity_score = self.substitution_matrix.ix[str(self.seq_B.seq[row_index - 1]), str(self.seq_A.seq[column_index - 1])]

                    # Decide whether to use gap_extension or gap_opening penalty score
                    if self.extending_gap == False:
                        gap_penalty = self.gap_opening_penalty
                    else:
                        gap_penalty = self.gap_extension_penalty

                    # Calculate possible scores for current cell in scoring matrix
                    score_up = self.scoring_matrix[row_index - 1, column_index] + gap_penalty
                    score_left = self.scoring_matrix[row_index, column_index] + gap_penalty
                    score_diag = self.scoring_matrix[row_index - 1, column_index - 1] + similarity_score

                    # Take highest final score and assign to current in scoring matrix
                    final_score = max(score_up, score_left, score_diag)

                    # Toggle gap_opening of gap_extension
                    if max(score_up, score_left, score_diag) == score_diag:
                        self.extending_gap = False
                    else:
                        self.extending_gap = True

                    self.scoring_matrix[row_index, column_index] = final_score if final_score > 0 else 0

        assert len(np.argwhere(self.scoring_matrix == -1)) == 0, "The scoring matrix was not completely filled in!"
        # np.savetxt("SCORE_current.csv", self.scoring_matrix, delimiter=",")

    def traceback(self):
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

            print("NORMALIZE == TRUE")
            print("SHORTER LENGTH: {}".format(shorter_length))
            print("NORMALIZED SCORE: {}".format(current_score / shorter_length))

        else:
            print("NORMALIZE == FALSE")

            self.max_alignment_score.append(current_score)

        print("Current Score: {}".format(current_score))

        seq_A_alignment = []
        seq_B_alignment = []

        current_row_index = row_index_high
        current_column_index = column_index_high

        # Begin traceback, stop when (above | left | diagonal) == 0
        while self.scoring_matrix[current_row_index, current_column_index] > 0:

            # Get [column, row] index for cell (above | left | diagonal) to current cell with highest score
            left_score = self.scoring_matrix[current_row_index, current_column_index - 1]
            above_score = self.scoring_matrix[current_row_index - 1, current_column_index]
            diag_score = self.scoring_matrix[current_row_index - 1, current_column_index - 1]

            next_move = max(diag_score, left_score, above_score)
            # seq_A - represented across COLUMNS
            # seq_B - represented down ROWS

            if next_move == diag_score:
                seq_A_alignment.append(str(self.seq_A.seq[current_column_index - 1]))
                seq_B_alignment.append(str(self.seq_B.seq[current_row_index - 1]))
                current_row_index -= 1
                current_column_index -= 1

            elif next_move == left_score:
                seq_A_alignment.append(str(self.seq_A.seq[current_column_index - 1]))
                seq_B_alignment.append('-')
                current_column_index -= 1

            elif next_move == above_score:
                seq_A_alignment.append('-')
                seq_B_alignment.append(str(self.seq_B.seq[current_row_index - 1]))
                current_row_index -= 1

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