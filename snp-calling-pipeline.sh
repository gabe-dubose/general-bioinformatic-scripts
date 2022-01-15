#!/bin/bash

<<COMMENT
This script implements an all-in-one pipeline that identifies single nucleotied polymorphisms and
indel mutations in given FASTQ files. BWA-mem is used for alighment, GATK is used to improve alignment, and
samtools is used for further processing and to call variants. Final output is a VCF file.
COMMENT

realign=0
gunzip=0
v=0
index=0
help=0

while getopts "a:b:r:e:o:f:zvih" option
do
case $option in
        a) reads1=$OPTARG;;
        b) reads2=$OPTARG;;
        r) ref=$OPTARG;;
        e) realign=1;;
        o) output=$OPTARG;;
        f) millsfile=$OPTARG;;
        z) gunzip=1;;
        v) v=1;;
        i) index=1;;
        h) help=1;;
    esac
done


#Echo help option
if [ $help -eq 1 ]; then
    echo "VC pipeline usage:
    -a reads1 file in FASTQ format (required)
    -b reads2 file in FASTQ format (required)
    -r reference sequence in FASTA format (required)
    -e realign your data to reduce the number of miscalled INDELS (optional, default=NULL)
    -o output identified: tag used to identify output files, will go before file extension (required)
    -f Mills file location (required)
    -z output VCF file is gzipped (*.vcf.gz) (optional, default=NULL)
    -i index output BAM file (optional, default=NULL)
    -h help - see command line options (optional, default=NULL)"
    exit 1
fi

#Check to see if input files exist
error=$(find $ref 2>&1)
error2=$(find $reads1 2>&1)
error3=$(find $reads2 2>&1)
error4=$(find $millsfile 2>&1)
if [ "$error" != "$ref" ]; then
    echo $ref was not found
    exit 1
elif [ "$error2" != "$reads1" ]; then
    echo $file2 was not found
    exit 1
elif [ "$error3" != "$reads2" ]; then
    echo $reads2 was not found
    exit 1
elif [ "$error4" != "$millsfile" ]; then
    echo $millsfile was not found
    exit 1
else
    echo "All input files have been located"
fi

#See if output vcf file exists
error4=$(find ${output}.vcf.gz 2>&1)
if [ "$error4" == "${output}.vcf.gz" ]; then
    echo "VCF file already exists. Would you like to write over it?(y/n)"
    read overwrite
    if [ "$overwrite" == "y" ]; then
        echo "${output}.vcf.gz will be overwritten."
        rm ${output}.vcf.gz
    elif [ "$overwrite" == "n" ]; then
        echo "Process will be stopped."
        exit 1
    fi
fi


#Initial alignment 
#index reference sequence
if [ $v -eq 1 ]; then
    echo "Indexing reference sequence..."
fi
bwa index $ref

#align reads to reference
if [ $v -eq 1 ]; then
    echo "Aligning reads to reference..."
fi
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > ${output}.sam

#clean read-pairing and flag information
if [ $v -eq 1 ]; then
    echo "Cleaning read-pairing and flagging information..."
fi
samtools fixmate -O bam ${output}.sam ${output}_fixmate.bam

#Convert cleaned .bam files to .sam files for input into sorting function
if [ $v -eq 1 ]; then
    echo "Converting BAM file to SAM file for sorting..."
fi
samtools view -h -o ${output}_fixmate.sam ${output}_fixmate.bam

#Sorts into coordinate order
if [ $v -eq 1 ]; then
    echo "Sorting into coordinate order..."
fi
samtools sort -O bam -o ${output}_sorted.bam -T /tmp/${output}_temp ${output}_fixmate.sam


#Realignment and improvement 
if [ $v -eq 1 ]; then
    echo "Beginning realignment..."
fi

#Create dictionary file from reference'
if [ $v -eq 1 ]; then
    echo "Creating dictionary file from reference..."
fi
java -jar picard.jar CreateSequenceDictionary R=$ref O=${ref}.dict
#rename the dictionary file because an error occurs when it is .fa.dict, as opposed to .dict
ls *.fa.dict | sed -e 'p;s/.fa.dict$/.dict/' | xargs -n2 mv

#create index file for reference
if [ $v -eq 1 ]; then
    echo "Indexing reference..."
fi
samtools faidx $ref

#create index file for sorted .bam file
if [ $v -eq 1 ]; then
    echo "Cleaning read-pairing and flagging information..."
fi
samtools index -b ${output}_sorted.bam

#Create target intervals file using sorted and indexed BAM file and VCF file of known indels
java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I ${output}_sorted.bam -o ${output}.intervals --known $millsfile -log ${output}_target_creator.log

#Realign indels using target intervals 
if [ $realign -eq 1 ]; then
    echo "Realigning indels..."
    java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I ${output}_sorted.bam -targetIntervals ${output}.intervals -known $millsfile -o ${output}_realigned.bam -log ${output}_indel_realignment.log
fi

#Format into 1 log file
cat ${output}_target_creator.log ${output}_indel_realignment.log >> ${output}.log
rm ${output}_target_creator.log ${output}_indel_realignment.log

#Index final BAM indel file (optional)
if [ $realign -eq 1 ]; then
    if [ $index -eq 1 ]; then
        echo "Indexing realigned BAM file..."
        samtools index ${output}_realigned.bam
    fi
elif [ $realign -eq 0 ]; then
    if [ $index -eq 1 ]; then
        echo "Indexing realigned BAM file"
        samtools index ${output}_sorted.bam
    fi
fi


#Variant call
if [ $v -eq 1 ]; then
    echo "Writing variant call file..."
fi

if [ $realign -eq 1 ]; then
    bcftools mpileup -Ou -f $ref ${output}_realigned.bam | bcftools call -vmO z -o ${output}.vcf.gz
elif [ $realign -eq 0 ]; then
    bcftools mpileup -Ou -f $ref ${output}_sorted.bam | bcftools call -vmO z -o ${output}.vcf.gz
fi

#Filter low quality calls
if [ $v -eq 1 ]; then
    echo "Filtering variant call file..."
fi
bcftools filter -O z -o ${output}_filtered.vcf.gz -s LOWQUAL -i'%QUAL>10' ${output}.vcf.gz

#Process vcf file into bed file
if [ $v -eq 1 ]; then
    echo "Convering VCF file into BED format..."
fi


gunzip ${output}_filtered.vcf.gz
#subtract alt from ref so that deletions are negative (makes more sense in my head)
sed '/^#/d' ${output}_filtered.vcf | cut -f1,2,4,5 | awk ' {print length($4) - length($3)} ' > len.txt
sed '/^#/d' ${output}_filtered.vcf | cut -f1 | sed 's/^...//' > nums.txt
sed '/^#/d' ${output}_filtered.vcf | cut -f2 > pos.txt
paste pos.txt len.txt | awk ' {print ($1) + ($2)} ' > alt.txt
paste nums.txt pos.txt alt.txt len.txt > variants.txt
rm pos.txt len.txt nums.txt alt.txt

#separate snps from indels
awk '{if ($4~/0/) {print $0}}' variants.txt > snps.txt
awk '{if ($4!~/0/) {print $0}}' variants.txt > indels.txt

if [ $gunzip -eq 1 ]; then
    gzip ${output}_filtered.vcf
fi