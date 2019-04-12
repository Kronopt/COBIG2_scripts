#!/bin/bash

# run fastqc on sample fastq.gz, on T number of threads

if ls | grep ".fastq.gz" -q; then  # looks for file names containing fastq.gz
	echo "How many Threads?"
	read T
	mkdir -p fastqc_OUT
	ls | grep "fastq.gz" | xargs -n$T fastqc -o fastqc_OUT -t $T
else
	echo "fastq.gz files not found"
fi
