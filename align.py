# CS-167 Global Sequence Alignment
# author: Jeffrey Lin
# Date: 2/7/2023
#
# Purpose: Determines the optimal global sequence alignment of two sequences 
#          within a FASTA file, according to the scoring scheme provided by the 
#          user 
# 

import numpy as np
import sys as sys

# parse_FASTA_file
# purpose: Finds and returns the two sequences within the FASTA file  
# input: FASTA file containing the sequences
# output: returns the two sequences as strings x and y respectively
# Effects: none

def parse_FASTA_file(file):
    Sequence_counter = 0
    x = ""
    y = ""

    #reads in sequences into strings
    for line in file:
        if ">" in line:
            Sequence_counter += 1
            continue
        elif Sequence_counter == 1:
            x += line.rstrip()
        elif Sequence_counter == 2:
            y += line.rstrip()

    return x, y 

# Purpose: The GlobalSeqAlignment returns the optimal global alignments of 
#          two sequences
class GlobalSeqAlignment:

    # Constructor
    # purpose: Initializes an instance of GlobalSeqAlignment class  
    # input: sequences x adn y, 3 integers values for scoring scheme 
    # output: none
    # Effects: creates instance of GlobalSeqAlignment Class
    def __init__(self, x, y, M, m, g):
        self.x = x
        self.y = y
        self.M = M
        self.m = m
        self.g = g

        # dimensions of the global sequence alignment matrix
        self.dimX = len(x) + 1
        self.dimY = len(y) + 1

        # create matrix that stores alignment scores  
        self.alignScore = np.empty([self.dimY, self.dimX], dtype = int)

        # stack to keep track of traceback for optimal global alignment 
        self.traceback_stack = []
        # matrix that has traceback for all values in the alignment matrix
        self.traceback_matrix = np.empty([self.dimY, self.dimX], dtype = tuple)

    # get_optimal_alignment
    # purpose: finds the optimal globe alignment of two sequences
    # input: none 
    # output: returns the optimal alignments for sequences x and y
    # Effects: none
    def get_optimal_alignment(self):

        self.__initialize()
        self.__fill_matrix()
        self.__fill_traceback_stack()

        return self.__get_traceback()

    # initialize
    # purpose: fills first row, column of scoring matrix based on gap penalty
    # input: none 
    # output: none
    # Effects: none
    def __initialize(self):

        #initializes first row of the matrix according to the gap penalty  
        for i in range(self.dimX):
            self.alignScore[0][i] = self.g * i
            self.traceback_matrix[0][i] = (0, i, "left")
        
        #initializes first column of the matrix according to the gap penalty
        for j in range(self.dimY):
            self.alignScore[j][0] = self.g * j
            self.traceback_matrix[j][0] = (j, 0, "up")
        
        self.traceback_matrix[0][0] = (0, 0,"start")

    # fill_matrix
    # purpose: fills in scoring matrix according to Needleman-Wunsch algorithm
    # input: none 
    # output: none
    # Effects: initializes the rest of the scoring matrix
    def __fill_matrix(self):

        #iterates through and fills entire matrix
        for i in range(1, self.dimX):
            for j in range(1, self.dimY):
                self.alignScore[j][i] = self.__score_alignments(i, j)
        
    # score_alignments
    # purpose: Determines which of 3 traceback directions yields best alignment
    # input: index x and index y of scoring matrix
    # output: none
    # Effects: none
    def __score_alignments(self, index_X, index_Y):

        seq_X_idx = index_X - 1
        seq_Y_idx = index_Y - 1

        # checks if substrings of sequences match
        if self.x[seq_X_idx] == self.y[seq_Y_idx]:
            Score = self.M
        else:
            Score = self.m

        # performs 3 calculations from Needleman-Wunsch algorithm
        diag = self.alignScore[index_Y - 1][index_X - 1] + Score
        left = self.alignScore[index_Y - 1][index_X] + self.g
        up = self.alignScore[index_Y][index_X - 1] + self.g

        # determines the max of the 3 scores
        max_score, traceback_pointer = self.__max_score(diag, left, up)

        # populates traceback matrix according to the index of the max score
        self.__populate_traceback(traceback_pointer, index_X, index_Y)

        return max_score
    
    # max_score
    # purpose: determines which of the 3 tracebacks has the greatest score
    # input: scores for tracing backwards diagonally, to the left, and upwards
    # output: returns the max of the 3 tracebacks and a string representing 
    #         the direction of the max traceback
    # Effects: none

    def __max_score(self, diag, left, up):

        # determines which of the 3 tracebacks has the max score and returns it
        if (diag == left and left == up):
            return diag, "diag"
        elif (diag == left and diag > up and left > up):
            return diag, "diag"
        elif (diag == up and diag > left and up > left):
            return diag, "diag"
        elif (diag < left and diag < up):
            traceback = self.__traceback_helper(diag, left,  up)
            return max(left, up), traceback
        elif(diag < left and diag < up and left == up):
            return up, "up"
        else:
            traceback = self.__traceback_helper(diag, left, up)
            return max(diag, left, up), traceback

    # traceback_helper
    # purpose: in cases where the max function is used, determines which value
    #          was the largest and returns string corresponding to that value
    # input: scores for tracing backwards diagonally, to the left, and upwards
    # output: returns a string representing 
    #         the direction of the max traceback
    # Effects: none

    def __traceback_helper(self, diag, left, up):
        
        # returns string corresponding to max value from max function
        if (diag < left and diag < up):
            if (left > up):
                return "left"
            else:
                return "up"
        else:
            if (diag > left and diag > up):
                return "diag"
            elif (left > diag and left > up):
                return "left"
            else:
                return "up"

    # populate_traceback
    # purpose: fills each index in traceback matrix with tuple showing where
    #          each index traces back to
    # input: string representing direction of traceback, indices x and y
    # output: none
    # Effects: initializes traceback matrix
    def __populate_traceback(self, traceback_pointer, idx_X, idx_Y):
        
        # fills traceback matrix with indices corresponding to diag, left, up
        if (traceback_pointer == "diag"):
            self.traceback_matrix[idx_Y][idx_X] = (idx_Y - 1, idx_X - 1, "diag")
        elif (traceback_pointer == "left"):
            self.traceback_matrix[idx_Y][idx_X] = (idx_Y -1, idx_X, "left")
        else:
            self.traceback_matrix[idx_Y][idx_X] = (idx_Y, idx_X - 1, "up")

    # fill_traceback_stack
    # purpose: wrapper function for recurisve stack_helper function
    # input: none
    # output: none
    # Effects: none
    def __fill_traceback_stack(self):
        len_X = len(self.x)
        len_Y = len(self.y)
        
        self.__stack_helper(len_X, len_Y)

    # stack_helper
    # purpose: fills stack with the tracebacks of the optimal global alignment
    # input: the length of sequences x and y
    # output: none
    # Effects: initializes traceback_stack

    def __stack_helper(self, idx_X, idx_Y):
            
        # fills stack with traceback until reaching index (0,0)
        if (self.traceback_matrix[idx_Y][idx_X][2] == "start"):
            return 
        else:
            prev = self.traceback_matrix[idx_Y][idx_X]
            self.traceback_stack.insert(0, prev)

            return self.__stack_helper(prev[1], prev[0])

    # get_traceback
    # purpose: forms optimal alignments of two sequences
    # input: none
    # output: returns tuple with optimal alignments of sequences x and y
    # Effects: creates tuple with optimal alignments of sequence x and y

    def __get_traceback(self):

        l2 = []
        sequence_X = ""
        sequence_Y = ""

        # continues to pop tracebacks from stack until it is empty    
        while not self.traceback_stack == l2:
            temp_tuple = self.traceback_stack[0]
            self.traceback_stack.pop(0)
    
            # append the sequences with optimal alignments of substrings
            if temp_tuple[2] == "diag":
                sequence_X += self.x[temp_tuple[1]]
                sequence_Y += self.y[temp_tuple[0]]
            elif temp_tuple[2] == "up":
                sequence_X += self.x[temp_tuple[1]]
                sequence_Y += "-"
            elif temp_tuple[2] == "left":
                sequence_X += "-"
                sequence_Y += self.y[temp_tuple[0]]
            elif temp_tuple[2] == "start":
                break
                
        return (sequence_X, sequence_Y)

    
def main():
    x,y = parse_FASTA_file(sys.stdin)
    
    # determine scoring scheme
    if len(sys.argv) == 4: 
        Match = int(sys.argv[1])
        Mismatch = int(sys.argv[2])
        Gap = int(sys.argv[3])
    elif len(sys.argv) == 1:
        Match = 4
        Mismatch = -2
        Gap = - 2

    align = GlobalSeqAlignment(x, y, Match , Mismatch, Gap)
    x_prime, y_prime = align.get_optimal_alignment()

    print(x_prime)
    print(y_prime)

if __name__ == "__main__":
    main()


