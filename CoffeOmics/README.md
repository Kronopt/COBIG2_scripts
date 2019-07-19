## Scripts related to the analysis of a RNAseq dataset

|  |  |
|---------|--------------|
| **all_blast2go.bash** | Runs Blast2GO for every FASTA file in a folder (recursively) |
| **bai_create.bash** | Creates bai index files for each mapped .bam file
| **build_fasta_from_gene_names.py** | Builds 2 FASTA files (up-regulated and down-regulated genes) with sequences that correspond to genes present in the output of DESEQ2 analysis |
| **descendant_goterms.py** | Identifies descendant GO terms and verify which genes are annotated with those descendants |
| **DESeq2_EdgeR_multiple** | Performs Differential Expression analysis with DESeq2 and EdgeR |
| **DESeq2_multiple.r** | Performs Differential Expression analysis with DESeq2 |
| **fastqc_multiple.bash** | Runs FastQC on sample 'fastq.gz', on 'T' number of threads |
| **gff_extract_all_swissprotids.py** | Extracts SwissProt IDs of each gene from a gff3 functional annotation file |
| **gmt_generate.py** | Generates .gmt and enrichments files suitable to be used in Cytoscape's EnrichmentMap |
| **htseq-count_multiple.bash** | htseq-counts multiple assemblies |
| **mean_logfc.py** | Generates a table of the mean LogFC for each DEG results file |
| **merge_fastq.gz.bash** | Merges sample file parts into single sample fastq.gz files |
| **PCA_htseq.R** | Performs PCA analysis of htseq-count data |
| **ranked_list_for_enrichment_map.py** | Generates ranked lists of genes (gene_name, FDR_padj) for use in a GSEA |
| **samtools_stats_multiple.bash** | samtools flagstat/idxstats/stats and plot-bamstats multiple assemblies |
| **STAR_multiple.bash** | Maps multiple RNA-seq samples with STAR |
| **top_10_DEGs.py** | Generates tables of the top 10 up- and down-regulated DEGs for each comparison |
| **top_10_goterms.py** | Generates a table of top 10 GO-terms |
| **trimmomatic_multiple.bash** | Trims multiple fastqc.gz files |
| **UniprotIdRetrieval.py** | Downloads protein sequences from Uniprot (adapted from UniprotIdRetrieval) |
| **venn.py** | Generates venn diagrams for total expressed, up- and down-regulated genes |
