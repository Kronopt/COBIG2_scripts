#!/bin/bash

# samtools flagstat/idxstats/stats and plot-bamstats multiple assemblies

if ls | grep "fastq.gz_RUN" -q; then

	echo "How many Threads?"
	read T

	# for each assembly folder
	for assembly_folder in $(ls | grep "fastq.gz_RUN"); do
		mkdir -p $assembly_folder/stats_counts/stats_plots
		bam_file=$(ls $assembly_folder | grep ".bam$")

		samtools flagstat -@ $T $assembly_folder/$bam_file \
			> $assembly_folder/stats_counts/$bam_file.flagstat

		samtools idxstats -@ $T $assembly_folder/$bam_file \
			> $assembly_folder/stats_counts/$bam_file.idxstats

		samtools stats -@ $T $assembly_folder/$bam_file \
			> $assembly_folder/stats_counts/$bam_file.stats

		plot-bamstats -p $assembly_folder/stats_counts/stats_plots/ \
			$assembly_folder/stats_counts/$bam_file.stats
	done
else
	echo "Assembly folders not found"
fi
