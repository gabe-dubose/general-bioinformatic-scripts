#!/usr/bin/env python3

'''
This script uses BLAST to identify orthologous genes between two input files in FASTA format.
Command line BLAST functions are wrapped for alignment, and the user may specify -p for a protein file
or -n for a nucleotide file, which will tell the program to use BLASTp or BLASTn, respectively.
'''

import argparse
import subprocess
import shlex

#Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i1', '--file1', required = False,
    metavar = "Sequence1", help="FASTA nucleotide or protein file")
parser.add_argument('-i2', '--file2', required = False,
    metavar = "Sequence2", help="FASTA nucleotide or protein file")
parser.add_argument('-t', '--sequence_type', required = False,
    metavar = "Sequence type", help="Either nucleotide(n) or protein(p)")
parser.add_argument('-o', '--output_file', required = False,
    metavar = "Output file name", help="The name of the output file")
args = parser.parse_args()

#determine which blast type to use
dbtype = ''
blast = ''
if args.sequence_type == 'n':
    dbtype = 'nucl'
    blast = 'blastn'
elif args.sequence_type == 'p':
    dbtype = 'prot'
    blast = 'blastp'


#create blast database for each input file
subprocess.call(shlex.split(f"makeblastdb -in {args.file1} -dbtype {dbtype} -out db1"))
subprocess.call(shlex.split(f"makeblastdb -in {args.file2} -dbtype {dbtype} -out db2"))


#blast query sequences
subprocess.call(shlex.split(f"{blast} -query {args.file1} -db db2 -outfmt 6 -out 1v2.txt"))
subprocess.call(shlex.split(f"{blast} -query {args.file2} -db db1 -outfmt 6 -out 2v1.txt"))

#Remove blast database file: for some reason rm db1* db2* didn't work?
#subprocess.check_output("rm db.*, shell=True") <- will get the job done
subprocess.call(shlex.split(f"rm db1.nhr db1.nin db1.nsq db2.nhr db2.nin db2.nsq"))

#open and read blast output files
file1 = open('1v2.txt', 'r')
file2 = open('2v1.txt', 'r')
line1 = file1.readlines()
line1.sort()
line2 = file2.readlines()
line2.sort()

#get columns of interest and store in lists for blast1v2
query1 = []
evalues1 = []
hit1 = []
for column in line1:
    e = float((column.split('\t')[10]))
    evalues1.append(float(e))
    q = column.split('\t')[0]
    query1.append(q)
    b = column.split('\t')[1]
    hit1.append(b)

#get columns of interest and store in lists for blast2v1
query2 = []
evalues2 = []
hit2 = []
for column in line2:
    e = float((column.split('\t')[10]))
    evalues2.append(float(e))
    q = column.split('\t')[0]
    query2.append(q)
    b = column.split('\t')[1]
    hit2.append(b)


#find best hit, store in dictionary with structure dict = {'query sequence': 'target sequence'} for blast1v2
best_dict1 = {}
for i in range(len(query1)-1):
    if query1[i] == query1[i+1] and evalues1[i] > evalues1[i+1]:
        best_dict1[query1[i]] = str(hit1[i])
    elif query1[i] == query1[i+1] and evalues1[i] < evalues1[i+1]:
        best_dict1[query1[i]] = str(hit1[i+1])
    elif query1[i] != query1[i+1]:
        best_dict1[query1[i]] = str(hit1[i])

#find best hit, store in dictionary with structure dict = {'query sequence': 'target sequence'} for blast2v1
best_dict2 = {}
for i in range(len(query2)-1):
    if query2[i] == query2[i+1] and evalues2[i] > evalues2[i+1]:
        best_dict2[query2[i]] = str(hit2[i])
    elif query2[i] == query2[i+1] and evalues2[i] < evalues2[i+1]:
        best_dict2[query2[i]] = str(hit2[i+1])
    elif query2[i] != query2[i+1]:
        best_dict2[query2[i]] = str(hit2[i])

#iterate through best hits 1 and 2 dictionaries and write orthologs to outfile
ortholog_count = 0
for hit1 in best_dict1:
    for hit2 in best_dict2:
        if best_dict1[hit1] == hit2 and best_dict2[hit2] == hit1:
            ortholog_count += 1
            with open(f'{args.output_file}_find_ortholog.output', 'a') as orthologs:
                orthologs.write(f'{hit1}\t{best_dict1[hit1]}\n')

#make summary file
with open(f'{args.output_file}_README.txt', 'a') as readme:
    readme.write(f'Input 1 BLAST hits against input 2: {len(query1)}\n')
    readme.write(f'\tHits after filtering for best hits only: {len(best_dict1)}\n')
    readme.write(f'Input 2 BLAST hits against input 1: {len(query2)}\n')
    readme.write(f'\tHits after filtering for best hits only: {len(best_dict2)}\n')
    readme.write(f'Number of orthologous genes: {ortholog_count}\n')

#delete initial blast output files
subprocess.call(shlex.split(f"rm 1v2.txt 2v1.txt"))