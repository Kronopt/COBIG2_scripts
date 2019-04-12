## Scripts related to the analysis of a RNAseq dataset

|  |  |
|---------|--------------|
| **all_blast2go.bash** | Runs Blast2GO for every FASTA file in a folder (recursively) | 
| **bai_create.bash** | Creates bai index files for each mapped .bam file
| **build_fasta_from_gene_names.py** | Builds 2 FASTA files (up-regulated and down-regulated genes) with sequences that correspond to genes present in the output of DESEQ2 analysis |
| **DESeq2_multiple.r** | Performs Differential Expression analysis with DESeq2 |
| **fastqc_multiple.bash** | Runs FastQC on sample 'fastq.gz', on 'T' number of threads |
| **gff_extract_all_swissprotids.py** | Extracts SwissProt IDs of each gene from a gff3 functional annotation file |
| **htseq-count_multiple.bash** | htseq-count multiple assemblies |
| **merge_fastq.gz.bash** | Merge sample file parts into single sample fastq.gz files |
| **PCA_htseq.R** | Performs PCA analysis of htseq-count data |
| **s1_table_generate.py** | Generates table of top 10 GO-terms |
| **samtools_stats_multiple.bash** | samtools flagstat/idxstats/stats and plot-bamstats multiple assemblies |
| **STAR_multiple.bash** | Map multiple RNA-seq samples with STAR |
| **top_10_DEGs.py** | Generates tables of the top 10 up- and down-regulated genes for each comparison |
| **trimmomatic_multiple.bash** | Trim multiple fastqc.gz files |
| **UniprotIdRetrieval.py** | Downloads protein sequences from Uniprot (adapted from UniprotIdRetrieval) |
| **venn.py** | Generates venn diagrams for total expressed, up- and down-regulated genes |
