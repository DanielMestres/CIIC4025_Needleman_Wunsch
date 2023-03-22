# Dynamic Programming Needleman-Wunsch Algorithm
# Developed following approach in Biological Sequence Alignments Slides 

import sys
import csv

# Scoring Matrix, match = +1, mismatch  = -1, d = gap = -2
def Scoring_Matrix(A, B):
    if(A == B):
        return match
    else:
        if(A != B):
            return mismatch
        else: 
            return d

# Needleman-Wunsch implementation
def NW_algo(seq1, seq2):
    # Initialize a m*n matrix
    m = len(seq1)
    n = len(seq2)
    matrix = [[0 for y in range(n+1)] for x in range(m+1)]

    # Init Step
    for j in range(n+1):
        matrix[0][j] = d * j
    for i in range(m+1):
        matrix[i][0] = d * i
    
    # Get max scores and fill matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            matrix[i][j] = max(matrix[i-1][j-1] + Scoring_Matrix(seq1[i-1], seq2[j-1]), matrix[i][j-1] + d, matrix[i-1][j] +d)

    # Backtracking step
    res_eq1 = ""
    res_eq2 = ""
    while(m > 0 or n > 0):
        if(matrix[m][n] == (matrix[m-1][n] + d)):
           res_eq1 = seq1[m-1] + res_eq1
           res_eq2 = "-" + res_eq2
           m -= 1
        else:
            if(matrix[m][n] == (matrix[m][n-1] + d)):
                res_eq1 = "-" + res_eq1
                res_eq2 = seq2[n-1] + res_eq2
                n -= 1
            else: 
                res_eq1 = seq1[m-1] + res_eq1
                res_eq2 = seq2[n-1] + res_eq2
                m -= 1
                n -= 1

    return(res_eq1, res_eq2)
            

# Main
match = 1
mismatch = -1
d = -2

if(len(sys.argv) > 1):
    with open(sys.argv[1]) as csv_file:
              reader = csv.DictReader(csv_file)
              for row in reader:
                  seq1 = row['sequence1']
                  seq2 = row['sequence2']
                  res1, res2 = NW_algo(seq1, seq2)
                  print(res1, res2)
else:
    print("Please provide a CSV input file")
