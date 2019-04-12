#!/bin/bash

# Trim multiple fastqc.gz files

if ls | grep ".fastq.gz" -q; then
	mkdir -p fastq_trimmed
	echo "How many Threads?"
	read T
	for input_file in $(ls | grep ".fastq.gz")
	do
		trimmomatic SE -threads $T $input_file fastq_trimmed/trimmed_$input_file ILLUMINACLIP:/home/pedrod/miniconda3/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:12 SLIDINGWINDOW:4:15 MINLEN:38
	done
else
	echo "files not found"
fi
