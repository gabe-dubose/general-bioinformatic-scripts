#!/usr/bin/env python3

'''
This script implements the Needleman-Wunsch alogrithm for sequence alignment. 
This was more for practice than anything.
'''

import sys

seq1 = sys.argv[1]
seq2 = sys.argv[2]

#initialize penalties and rewards
gap = -1
mis = -1
match = 1
#read each sequence file, store in separate lists

#open file, read second line
seq1 = open(f'{seq1}', 'r')
seq1 = seq1.readlines()
seq = seq1[1].rstrip()

#initialize list to store sequence and add every base to that list
sequence1 = []
for chr in seq:
    sequence1.append(chr)

#do the same as above for second sequence
seq2 = open(f'{seq2}', 'r')
seq2 = seq2.readlines()
seq = seq2[1].rstrip()
sequence2 = []
for chr in seq:
    sequence2.append(chr)

#set up scoring matrix,
#make first row in matrix and fill with gap penalties, 
# fill start of every row with gap penalty
matrix = []
for i in range(len(sequence2)+1):
    matrix.append([])
for i in range(len(sequence1)+1):
    matrix[0].append(gap*i)
x = 1
while x < len(matrix):
    matrix[x].append(x*gap)
    x+=1

#set up traceback matrix same as before
traceback = []
for i in range(len(sequence2)+1):
    traceback.append([])
for i in range(len(sequence1)+1):
    traceback[0].append(gap*i)
x = 1
while x < len(traceback):
    traceback[x].append(x*gap)
    x+=1

#matrix filling
for i in range(len(sequence2)):
    for j in range(len(sequence1)):
        diagonal = matrix[i][j]
        left = matrix[i+1][j] + gap
        up = matrix[i][j+1] + gap
        #compare adjacent values
        #see if diagonal is hit or miss
        if sequence2[i] == sequence1[j]:
            diagonal = diagonal + match
        elif sequence2[i] != sequence1[j]:
            diagonal = diagonal + mis
        #compare scores, add highest to matrix
        if diagonal >= up and diagonal >= left:
            matrix[i+1].append(diagonal)
            traceback[i+1].append('d')
        elif up >= left and up > diagonal:
            matrix[i+1].append(up)
            traceback[i+1].append('u')
        elif left > up:
            matrix[i+1].append(left)
            traceback[i+1].append('l')

#backtracking
#get length of longest sequence for iteration number
if len(sequence1) > len(sequence2):
    length = len(sequence1)
elif len(sequence2) > len(sequence1):
    length = len(sequence2)

#define lists that will hold alignment information
line1 = []
line2 = []
line3 = []

#loop through traceback matrix and find path of best alignment, 
#add bases and gap/mismatch information to lists
x=-1
y=-1
for i in range(length):
    target = traceback[x][y]
    if target == 'd':
        line1.append(sequence1[y])
        if sequence1[y] == sequence2[x]:
            line2.append('|')
        else:
            line2.append('*')
        line3.append(sequence2[x])
        x-=1
        y-=1
    elif target == 'l':
        line1.append(sequence1[y])
        line2.append(' ')
        line3.append('-')
        x = x
        y -= 1
    elif target == 'u':
        line1.append('-')
        line2.append(' ')
        line3.append(sequence2[x])
        x -= 1
        y = y


#print results
str = ""
print(str.join(line1[::-1]))
print(str.join(line2[::-1]))
print(str.join(line3[::-1]))
print(f'Alignment Score: {matrix[-1][-1]}')
