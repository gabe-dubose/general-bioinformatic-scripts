#!/usr/bin/env python3

'''
This script wraps MUMmer3's dnadiff function to calculate pairwise average nucleotide identity (ANI)
between all input files and returns a matrix of each comparison. Users may specify the number of 
parallel processes that this program uses to increase runtime speed. 
'''
import argparse
import subprocess
import shlex
import os
import csv
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--output', metavar = '<Output file name>')
parser.add_argument('-t', '--threads', type=int, metavar = '<Number of processes>')
parser.add_argument('files', nargs='+', type=str)
args = parser.parse_args()

#read in sequences input and store in list
sample_matrix = []
for sequence in args.files:
    sample_matrix.append(sequence)
sample_matrix.sort()

#generate list of individual comparisons to make
cmp_list = []
x=0
y=1
for i in range(1000):
    if y < len(sample_matrix):
        cmp1 = sample_matrix[x]
        cmp2 = sample_matrix[y]
        cmp_list.append([])
        cmp_list[i].append(cmp1)
        cmp_list[i].append(cmp2)
        y+=1
    if y > len(sample_matrix)-1:
        x=0
        y=1
        sample_matrix = sample_matrix[1:]

#define function to run dnadiff with input from cmp_list 
def run_dnadiff(comparisons):
    cmp1 = comparisons[0]
    cmp2 = comparisons[1]
    subprocess.call(shlex.split(f"dnadiff -p {cmp1}_to_{cmp2} {cmp1} {cmp2}"))
    with open(f"{cmp1}_to_{cmp2}.report", 'r') as output:
        lines = output.readlines()
        ani = str(lines[18:19])
        ani = " ".join(ani.split())
        ani = ani.split(' ')[1]
    with open('temp_out.txt', 'a') as temp:
        temp.write(f"{cmp1}\t{cmp2}\t{ani}\n")

pool = Pool(args.threads)
pool.map(run_dnadiff, cmp_list)

os.system("rm *_to_*")



#write output
#make empty matrix full of 0s
sample_matrix = []
for sequence in args.files:
    sample_matrix.append(sequence)
sample_matrix.sort()

outmatrix = []
for i in range(len(sample_matrix)+1):
    outmatrix.append([])
outmatrix[0].append('')
for sample in sample_matrix:
    outmatrix[0].append(sample)
count=1
for sample in sample_matrix:
    outmatrix[count].append(sample)
    count+=1
count=1
for line in outmatrix[1:]:
    for i in range(len(sample_matrix)):
        line.append(0)

#add 100s for identical samples
count=1
for i in range(len(sample_matrix)):
    outmatrix[count][count] += 100
    count+=1

lines = []
with open('temp_out.txt', 'r') as temp:
    for line in temp:
        lines.append(line)
lines.sort()

for result in lines:
    ref = result.split('\t')[0]
    seq = result.split('\t')[1]
    ani = float(result.rstrip('\n').split('\t')[2])
    print(f"{ref}\t{seq}\t{ani}")
    x=1
    y=1
    for i in range(len(outmatrix)):
        for j in range(len(outmatrix[0])):
            try:
                if ref == outmatrix[0][x] and seq == outmatrix[y][0]:
                    outmatrix[y][x] += ani
                if x < len(outmatrix[0]):
                    x+=1
                if x == len(outmatrix[0]):
                    x=1
                    y+=1
            except:
                continue

os.system("rm temp_out.txt")

with open('temp_out.txt', 'wt') as output:
    tsv_writer = csv.writer(output, delimiter='\t')
    tsv_writer.writerows(outmatrix)

os.system(f"sed -E 's/100.*//g' temp_out.txt > {args.output}")

os.system("rm temp_out.txt")

