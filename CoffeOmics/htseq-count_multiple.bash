#!/bin/bash

# htseq-count multiple assemblies

if ls | grep "fastq.gz_RUN" -q; then
	echo "How many Threads?"
        read T

	# create stats_counts folders inside each assembly folder
	for assembly_folder in $(ls | grep "fastq.gz_RUN"); do
		mkdir -p $assembly_folder/stats_counts
	done

	# do parallel htseq-count for each bam file on each assembly folder
	find . -type f -name \*.bam | sort | parallel --no-notice --bar -P $T \
		'htseq-count -i gene_name -q -s reverse -f bam {} \
		/home/pedrod/COFFEE/data/c_canephora_genome/c_canephora_gene_structural_and_functional_annotation.gtf \
		> {//}/stats_counts/{/}.htseq_counts'
else
	echo "Assembly folders not found"
fi
