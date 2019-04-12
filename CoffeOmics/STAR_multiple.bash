#!/bin/bash

# STAR map multiple samples

# SCRIPT PARAMETERS
# $1 = Samples FastQ's folder FULL PATH
# $2 = STAR Genome indexes folder FULL PATH
# $3 = Threads

# Parameters must be supplied
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
	echo "3 arguments needed:"
	echo "1 - Path to sample fastq.gz's"
	echo "2 - Path to STAR genome indexes"
	echo "3 - threads"
	exit 1
else  # 3 arguments were supplied
	if [ -z "$(ls -A)" ]; then  # if current folder is empty
		# list samples
		fastq_list=$(ls $1 | grep "fastq.gz" | sed s,^,$1,)

		# load genome index to memory
		STAR --genomeLoad LoadAndExit --genomeDir $2

		# ASSEMBLY
		for fastq in $fastq_list; do
			echo "Sample file path:" $fastq
			fastq_name=$(basename $fastq)
			mkdir $fastq_name\_RUN

			STAR --runThreadN $3 --readFilesIn $fastq --readFilesCommand zcat \
			--genomeDir $2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
			--limitBAMsortRAM 50000000000 --outFileNamePrefix $fastq_name\_RUN/$fastq_name\_
		done

		# unload genome index from memory
		STAR --genomeLoad Remove --genomeDir $2

	else
		echo "Folder is not empty. Folder must be empty."
	fi
fi
