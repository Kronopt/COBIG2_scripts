#!/usr/bin/env python3
# coding: utf-8

"""
Generates ranked lists of genes (gene_name, FDR_padj).
Intended for use in GSEA analysis followed by enrichment map generation.
"""

import csv, argparse, os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Generates ranked lists of genes (gene_name, FDR_padj).')
	parser.add_argument('seqtable_folder', help='folder containing seqtable.tab blast2go output files.')
	parser.add_argument('de_results_folder', help='folder containing results.csv DE output files.')
	parser.add_argument('output_folder', help='output folder for tables (.rnk files)')
	parser = parser.parse_args()

	# Handle folders
	if (os.path.exists(parser.seqtable_folder) and
		os.path.exists(parser.de_results_folder) and
		os.path.exists(parser.output_folder)):

		# seqtables
		print('Parsing seqtable.tab files ... ', end='')
		seqtab_dict = {}  # {1vs3_DOWN: [gene_name, ...], ... }
		list_seqtables = [seqtable_file for seqtable_file in os.listdir(parser.seqtable_folder)
							if 'seqtable.tab' in seqtable_file]
		for seqtable_file in list_seqtables:
			file_name = seqtable_file.split('_')[0] + '_' + seqtable_file.split('_')[2].split('.')[0]
			seqtab_dict[file_name] = []
			with open(os.path.join(parser.seqtable_folder, seqtable_file)) as seqtab_file:
				next(seqtab_file)  # ignore header
				for line in seqtab_file:
					seqtab_dict[file_name].append(line.split()[0])
		print('OK')

		# de_results
		print('Parsing DE results.csv files ... ', end='')
		list_de_results = [de_results_file for de_results_file in os.listdir(parser.de_results_folder)
							if 'results.csv' in de_results_file]
		de_results_dict = {}  # {1vs3: [(gene_name, FDR), ...], ...}
		for de_results_file in list_de_results:
			file_name = de_results_file.split('_')[0]
			de_results_dict[file_name] = []
			with open(os.path.join(parser.de_results_folder, de_results_file), newline='') as de_file:
				de_file_reader = csv.reader(de_file)
				next(de_file_reader)  # ignore header
				for row in de_file_reader:
					if row[6] != 'NA' and float(row[6]) < 0.1:  # FDR p-adjusted
						de_results_dict[file_name].append((row[0], row[6]))
		print('OK')

		# filter genes
		print('Filtering for genes in seqtables ... ', end='')
		final_genes = {}  # {1vs3: [(gene_name, FDR), ...], ...}
		for comparison in seqtab_dict:
			final_genes[comparison] = []
			for gene_name in seqtab_dict[comparison]:
				for other_gene_name, FDR in de_results_dict[comparison.split('_')[0]]:
					if gene_name == other_gene_name:
						final_genes[comparison].append((gene_name, FDR))
		print('OK')

		# write output file
		print('Writing output file ... ', end='')
		header = ['#', 'Gene name', 'FDR padj']
		for ranked_list in final_genes:
			ranked_list_name = ranked_list.split('_')
			with open(os.path.join(parser.output_folder,
				ranked_list_name[0] + '_ranked_list_' + ranked_list_name[1] + '.rnk'), 'w', newline= '') as csv_file:
				csv_writer = csv.writer(csv_file, delimiter='\t')
				csv_writer.writerow(header)  # write header
				for row in final_genes[ranked_list]:
					csv_writer.writerow(row)
		print('OK')

	else:
		print('Please input valid and existing paths for \'seqtable_folder\', \'de_results_folder\' and \'output_folder\'.')
