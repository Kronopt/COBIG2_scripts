#!/bin/bash

# merge sample file parts into single sample fastq.gz files
# 12 samples with 3 replicates each (A, B and C)

mkdir -p merged_fastqgz

for i in {1..12}
do
	for j in A B C
	do
		PATTERN="IICT_$i$j"  # looks for files with this name
		echo "Sample $i$j"
		if ls | grep $PATTERN -q; then
			ls | grep $PATTERN
			ls | grep $PATTERN | xargs cat > merged_fastqgz/$i$j.fastq.gz
			echo "done"
		else
			echo "files not found"
		fi
	done
done
