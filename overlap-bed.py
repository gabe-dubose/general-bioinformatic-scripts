#!/usr/bin/env python3

'''
This script takes two genomic coordinate files in BED format and identifies which elements
overlap wihtin a specified percent overlap threshold. The out put is a set of overlapping genomic
coordinates of the first file provided. However, the user can specify the --join option, whihc will return
the overlapping coordinates from both files.
'''

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i1', '--input_file1', metavar = "<Input file>", 
    help = "First input file in BED format")
parser.add_argument('-i2', '--input_file2', metavar = "<Input file>",
    help = "Second input file in BED format")
parser.add_argument('-m', '--minimum_overlap', metavar = "<Percentage>",
    help = "Minimum percent overlap", type = int)
parser.add_argument('-j', '--join', action = "store_true",
    help = "Output first set and second set of overlaps side by side")
parser.add_argument('-o', '--output', metavar = "<output_file_name>", 
    help = "Name of output file")
args = parser.parse_args()


def calculate_overlap(chrom, start1, stop1, start2, stop2):
    len = int(stop1) - int(start1)
    overlap = (min(int(stop1), int(stop2))) - (max(int(start1), int(start2)))
    if len > 0:
        percent_overlap = (overlap/len)*100
    if overlap > 0 and percent_overlap >= args.minimum_overlap:
        with open(args.output, 'a') as output:
            if args.join == True:
                output.write(f"{chrom}\t{start1}\t{stop1}\t{chrom}\t{start2}\t{stop2}\n")
            else:
                output.write(f"{chrom}\t{start1}\t{stop1}\n")
        print(f"{chrom}\t{start1}\t{stop1}")
     
dict1 = {}
with open(args.input_file1, 'r') as file1:
    for line in file1:
        if line.split('\t')[0] not in dict1.keys():
            dict1[line.split('\t')[0]] = []
            dict1[line.split('\t')[0]].append([line.split('\t')[1], line.rstrip().split('\t')[2]])
        else:
            dict1[line.split('\t')[0]].append([line.split('\t')[1], line.rstrip().split('\t')[2]])
file1.close()

dict2 = {}
with open(args.input_file2, 'r') as file2:
    for line in file2:
        if line.split('\t')[0] not in dict2.keys():
            dict2[line.split('\t')[0]] = []
            dict2[line.split('\t')[0]].append([line.split('\t')[1], line.rstrip().split('\t')[2]])
        else:
            dict2[line.split('\t')[0]].append([line.split('\t')[1], line.rstrip().split('\t')[2]])
file2.close()

for key in dict1:
    for i in range(len(dict1[key])):
        try:
            chrom = key
            start1 = dict1[key][i][0]
            stop1 = dict1[key][i][1]
            start2 = dict2[key][0][0]
            stop2 = dict2[key][0][1]
            if stop2 < start1:
                dict2[key] = dict2[key][1:]
            if stop1 < start2:
                continue
            for j in range(len(dict2)):
                start2 = dict2[key][j][0]
                stop2 = dict2[key][j][1]
                calculate_overlap(chrom, start1, stop1, start2, stop2)
        except:
            continue