## This code is borrowed  from somewhere in github.  I will add the credit once I find the author.

import sys
import math

#Sets the values for match, mismatch, and gap.
match = 1
mismatch = 0
gap = -1

def get_helper_data(seqA, seqB):
    seqB = seqB[:-1]
    seqA = "-"+seqA
    seqB = "-"+seqB

    match = 1
    mismatch = 0
    gap = -1

    row = len(seqA)
    col = len(seqB)
    return seqA, seqB, match, mismatch, gap, row, col

#Creates blank matrix
def create_matrix(row,col):
    A = [0] * row
    for i in range(row):
        A[i] = [0] * col
    return A

def isMatch(i,j, seqA, seqB):
    if seqA[i] == seqB[j]:
        matchVal = match
    else:
        matchVal = mismatch
    return matchVal

#Returns the new value if diagonal is used
def diag(A,i,j,seqA,seqB):
    return A[i-1][j-1] + isMatch(i,j, seqA, seqB)

#Returns the new value if up is used
def up(A,i,j,seqA,seqB):
    return A[i-1][j] + isMatch(i,j, seqA, seqB) + gap

#Returns the new value if left is used
def left(A,i,j,seqA,seqB):
    return A[i][j-1] + isMatch(i,j, seqA, seqB) + gap

#Fills matrix with correct scores.
def complete_matrix(A,row,col,seqA,seqB):
    for i in range(1,row):
        for j in range(1,col):
            A[i][j] = max(0,diag(A, i,j,seqA,seqB),up(A,i,j,seqA,seqB),left(A,i,j,seqA,seqB))
    return A

#FInd the highest scoring cell.
def get_max(A,row,col):
    local_max = [[0,0]]
    for i in range(row):
        for j in range(col):
            if A[i][j] == A[local_max[0][0]][local_max[0][1]]:
                local_max.append([i,j])
            elif A[i][j] > A[local_max[0][0]][local_max[0][1]]:
                local_max = [[i,j]]
    return local_max

#Gives you the next location.
def get_next(A,location):
    gap = -1
    i = location[0]
    j = location[1]
    maxVal = max(A[i-1][j-1],A[i-1][j]+gap,A[i][j-1]+gap)
    if A[i-1][j-1] == maxVal:
        return [i-1,j-1]
    #Is this the right ordering of the three?
    elif A[i][j-1]+gap == maxVal:
        return [i,j-1]
    else:
        return [i-1,j]

#Traces the path back given starting location
def trace_back(A,tracer):
    if A[tracer[len(tracer)-1][0]][tracer[len(tracer)-1][1]] == 0:
        return tracer
    next_cell = get_next(A,tracer[len(tracer)-1])
    #tracer.insert(0,next_cell)
    tracer.append(next_cell)
    return trace_back(A,tracer)

#Uses tracer to return final sequence
def get_seq(A,tracer,k,seq, seqA, seqB):
    if k == 0:
        original_sequence = seqA
    else:
        original_sequence = seqB
    N = len(tracer)
    for i in range(0,N-1):
        if tracer[i][k] == tracer[i+1][k]+1:
            seq = seq + original_sequence[tracer[i][k]]
        elif tracer[i][k] == tracer[i+1][k]:
            seq = seq + "-"
    return seq

#Shows the relevant lines for matching pairs
def get_middle(finalA,finalB):
    middle = ""
    for k in range(len(finalA)):
        mid = " "
        if finalA[k] == finalB[k]:
            mid = "|"
        middle = middle + mid
    return middle

def align_TE(seqA_backup, seqB_backup):
    seqA, seqB, match, mismatch, gap, row, col = get_helper_data(seqA_backup, seqB_backup)
    
    A = create_matrix(row,col)
    A = complete_matrix(A,row,col,seqA,seqB)

    num_answers = len(get_max(A,row,col))

    for i in range(num_answers):
        tracer = trace_back(A,[get_max(A,row,col)[i]])
        finalA = get_seq(A,tracer,0,"", seqA, seqB)
        finalB = get_seq(A,tracer,1,"", seqA, seqB)
        finalA = finalA[::-1]
        finalB = finalB[::-1]
        middle = get_middle(finalA,finalB)
        #print tracer
        #print("{}\n{}\n{}\n".format(finalA,middle,finalB))
        #print "finalA", finalA
        #print "finalB", finalB
        x = finalA.find('-')
        finalC = finalA[:x]
        break_point = seqA_backup.find(finalC)
        ##
        y = finalB.replace('-', '')
        z = seqB_backup.find(y)
        #print len(seqB_backup), len(finalB.replace('-', '')), y
        ##
        context = seqA_backup[:break_point - z]
        TE = seqA_backup[break_point - z:]
        return context, TE