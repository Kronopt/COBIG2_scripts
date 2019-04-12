#!/bin/bash

# create bai index files for each mapped bam

if ls | grep "_RUN" -q; then  # folders of mappings end in _RUN
	echo "How many Threads?"
	read T
	for mapping_folder in $(ls | grep "_RUN"); do
		cd $mapping_folder
		samtools index -@ $T $(ls | grep ".bam")
		cd ..
	done
else
	echo "Folders containing mapping outputs not found"
fi
