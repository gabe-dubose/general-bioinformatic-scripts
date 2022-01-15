#!/usr/bin/env python3

import argparse

'''
This script counts k-mers in a given protein or nucleotide sequence file in FASTA format. User can 
specify the k-mer lenght using the -k option. 
'''

#Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', required = False,
    metavar = "Sequence file", help="FASTA nucleotide or protein file")
parser.add_argument('-k', '--kmer_lenght', required = False,
    metavar = "k-mer lenght", help="The lenght of k-mers you want to count")
args = parser.parse_args()

seq_lines = []
lines = args.file.readlines()

print("Configuring file...")
#Get rid of first line
for number, line in enumerate(lines):
    if number not in [0]:
        seq_lines.append(line)

#format sequence into one big string
seq_lines = ''.join(seq_lines)
seq_lines = seq_lines.replace('\n', '')

#initialize dictionary that will hold each k-mer and its count
#initialize list that will store each k-mer in the sequence
kmer_dict={}
kmer_list = []

start = 0
end = int(kmer_length)

print(f"Creating {kmer_length}-mer dictionary...")
#get list of all kmers -> store in list
for i in range(len(seq_lines)):
    kmer = seq_lines[start:end]
    if len(kmer) < int(kmer_length):
        break
    kmer_list.append(kmer)
    start += 1
    end += 1
kmer_list.sort()

#configure dictionary: add every k-mer to dictionary, set value to 0, and erase duplicate values
#only adds unique k-mers because list is sorted (also I guess dictionaries only store unique values)
for kmer in kmer_list:
    kmer_dict[kmer] = 0

print(f"Counting {kmer_length}-mers...")
#count k-mers by comparing each dictionary element to the list elements
for key in kmer_dict:
    for kmer in kmer_list:
        if key == kmer:
            kmer_dict[kmer] += 1

for key in kmer_dict:
    print(f"{key}\t{kmer_dict[key]}")

print("Writing results...")
for key in kmer_dict:
    with open(f'{kmer_length}-mers.txt', 'a') as kmer_counts:
        kmer_counts.write(f"{key}\t{kmer_dict[key]}\n")