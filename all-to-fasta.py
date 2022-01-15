#!/usr/bin/env python3

'''
This script takes (and auto-detects) EMBL, FASTQ, GenBank, MEGA, SAM, and VCF files and converts them to a
FASTA file. The user may also specify the number of nucleotide characters to be included on each line. 
'''

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fold', required=False, default = 70, type = int,
    metavar = "<Line fold>", help = "The lenght of sequence per line; default=70")
parser.add_argument('-i', '--input', required=True,
    metavar = "<Input sequence file>", help = "Input sequence file in either EMBL, FASTQ, GenBank, MEGA, SAM, or VCF format")
args = parser.parse_args()

with open(args.input, 'r') as input_file:
    first_line = input_file.readline()
    second_line = input_file.readline()    
    #identify possible file types
    mega = bool(re.match("#MEGA", first_line))
    if mega == True:
        file_type = 'mega'
    fastq = bool(re.match("^@", first_line) and re.match("^[A-Za-z]", second_line))
    if fastq == True:
        file_type = 'fastq'
    embl = bool(re.match("^ID", first_line))
    if embl == True:
        file_type = 'embl'
    sam = bool(re.match("^@", first_line) and re.match("^@", second_line))
    if sam == True:
        file_type = 'sam'
    vcf = bool(re.match("^##", first_line))
    if vcf == True:
        file_type = 'vcf'
    genbank = bool(re.match("^LOCUS", first_line))
    if genbank == True:
        file_type = 'genbank'
    print(f'File type detected: {file_type}')
    
def wrap(string, length):
    return '\n'.join(string[i:i+length] for i in range(0, len(string), length) )

def convert_mega(file):
    sequence = ''
    with open(file, 'r') as input_file:
        for i in range(3):
            next(input_file)
        for line in input_file:
            if bool(re.match("^#", line)) == True:
                header = '>' + str(line).rstrip('\n').lstrip('#')
            else:
                sequence = sequence + str(line.rstrip('\n'))
        if bool(re.match("^[ATCGNatgcn]{4}", sequence)) == True:
            extension = 'fna'
        else:
            extension = 'faa'
    with open(f"{args.input}.{extension}",'a') as output:
        output.write(f'{header}\n{wrap(sequence, args.fold)}')

def convert_fastq(file):
    with open(file, 'r') as input_file:
        count = 0
        for line in input_file:
            count +=1
            if bool(re.match("^@", line)) == True:
                header = '>' + str(line).lstrip("@").rstrip('\n')
            if count % 2 == 0 and count % 4 != 0:
                sequence = line.rstrip('\n')
                if bool(re.match("^[ATCGNatgcn]{4}", sequence)) == True:
                    extension = 'fna'
                else:
                    extension = 'faa'
                with open(f"{args.input}.{extension}",'a') as output:
                    output.write(f'{header}\n{wrap(sequence, args.fold)}')

def convert_genbank(file):
    sequence = ''
    with open(file, 'r') as input_file:
        for line in input_file:
            if bool(re.match("^DEFINITION", line)) == True:
                header = '>' + line.lstrip("DEFINITION  ").rstrip("\n")
            if bool(re.search("[0-9]{1,} [A-Za-z]{10}", line)) == True:
                sequence = sequence + line.replace(' ','').lstrip('0123456789').rstrip('\n')
            if bool(re.match("^[ATCGNatgcn]{4}", sequence)) == True:
                extension = 'fna'
            else:
                extension = 'faa'
    with open(f"{args.input}.{extension}",'a') as output:
        output.write(f'{header}\n{wrap(sequence, args.fold)}')
       

def convert_embl(file):
    sequence = ''
    with open(file,'r') as input_file:
        for line in input_file:
            if bool(re.match("^ID", line)) == True:
                header = '>' + line.lstrip("ID   ").rstrip('\n')
            if bool(re.search("\s{5}[A-Za-z]{10}", line)) == True:
                sequence = sequence + line.replace(' ','').rstrip('0123456789\n ').lstrip('\s')
            if bool(re.match("^[ATCGNatgcn]{4}", sequence)) == True:
                extension = 'fna'
            else:
                extension = 'faa'
    with open(f"{args.input}.{extension}",'a') as output:
        output.write(f'{header}\n{wrap(sequence, args.fold)}')

def convert_sam(file):
    with open(file,'r') as input_file:
        for line in input_file:
            if bool(re.match("^@", line)) == True:
                next(input_file)
            else:
                header = '>' + line.split('\t')[0]
                sequence = line.split('\t')[9]
                if bool(re.match("^[ATCGNatgcn]{4}", sequence)) == True:
                    extension = 'fna'
                else:
                    extension = 'faa'
                with open(f"{args.input}.{extension}",'a') as output:
                    output.write(f'\n{header}\n{wrap(sequence, args.fold)}')
      
def convert_vcf(file):
    samples_dict = {}
    ref_seq = ''
    with open(file,'r') as input_file:
        for line in input_file:
            if bool(re.search("^#CHROM", line)) == True:
                samples = line.rstrip('\n').split('\t')[9:]
                for sample in samples:
                    samples_dict['>' + str(sample)] = []
            if bool(re.search("^#", line)) == False:
                ref_name = line.rstrip('\n').split('\t')[0]
                var = line.rstrip('\n').split('\t')[9:]
                ref = line.rstrip('\n').split('\t')[3]
                ref_seq += ref
                alt = line.rstrip('\n').split('\t')[4].split(',')
                for i in range(len(var)):
                    type = var[i][0]
                    sample = '>sample' + str(i+1)
                    if int(type) == 0:
                        samples_dict[sample].append(ref)
                    else:
                        samples_dict[sample].append(alt[int(type)-1])        
    samples_dict['>' + str(ref_name)] = ref_seq
    #print(samples_dict)
    for samples in samples_dict:
        header = samples
        sequence = ' '.join(samples_dict[samples])
        sequence = re.sub(' ', '', sequence)
        if bool(re.match("^[ATCGNatgcn]{4}", sequence)) == True:
            extension = 'fna'
        else:
            extension = 'faa'
        with open(f"{args.input}.{extension}",'a') as output:
            output.write(f'\n{header}\n{wrap(sequence, args.fold)}')
        #print(samples_dict[samples])

if file_type == 'mega':
    convert_mega(args.input)
if file_type == 'fastq':
    convert_fastq(args.input)
if file_type == 'genbank':
    convert_genbank(args.input)
if file_type == 'embl':
    convert_embl(args.input)
if file_type == 'sam':
    convert_sam(args.input)
if file_type == 'vcf':
    convert_vcf(args.input)
