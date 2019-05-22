#!/usr/bin/env python3
# coding: utf-8

"""
Generate .gmt and enrichments files suitable to be used in Cytoscape's EnrichmentMap.
"""

import argparse, os

def go_and_description(file_handle):
	final_dict = {}
	for line in file_handle:
		final_dict[line.split('\t')[0]] = line.rstrip()
	return final_dict

def enrichments(file_handle):
	header = next(file_handle)
	return header, go_and_description(file_handle)

def genes_folder(folder):
	final_dict = {}
	list_gene_files = [gene_file for gene_file in os.listdir(folder) if gene_file.startswith('GO_') and gene_file.endswith('.txt')]
	for gene_file in list_gene_files:
		with open(os.path.join(folder, gene_file)) as genes:
			go_term = gene_file.split('.')[0].replace('_', ':')
			if go_term not in final_dict:
				final_dict[go_term] = set()
			for gene in genes:
				final_dict[go_term].add(gene.rstrip())
	return final_dict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Generate .gmt and enrichments files.')
	parser.add_argument('enrichments_pos', type=argparse.FileType('r', encoding='UTF-8'),
						help='GSEA output .txt file for up-regulated genes. Should have the following columns:\n'
						'NAME, GS, GS DETAILS, SIZE, ES, NES, NOM p-val, FDR q-val, FWER p-val, RANK AT MAX, LEADING EDGE')
	parser.add_argument('enrichments_neg', type=argparse.FileType('r', encoding='UTF-8'),
						help='GSEA output .txt file for down-regulated genes. Should have the following columns:\n'
						'NAME, GS, GS DETAILS, SIZE, ES, NES, NOM p-val, FDR q-val, FWER p-val, RANK AT MAX, LEADING EDGE')
	parser.add_argument('go_and_description_up', type=argparse.FileType('r', encoding='UTF-8'),
						help='GO term + description file for up-regulated genes.\n'
						'File has no header and each line is in the form of <GO:####### \\t description>')
	parser.add_argument('go_and_description_down', type=argparse.FileType('r', encoding='UTF-8'),
						help='GO term + description file for down-regulated genes.\n'
						'File has no header and each line is in the form of <GO:####### \\t description>')
	parser.add_argument('genes_folder_up', help='folder containing up-regulated gene list files.\n'
						'Each file in this folder is named after a GO term (GO_#######) and each line of each file has a single gene name (Cc##_g#####)')
	parser.add_argument('genes_folder_down', help='folder containing down-regulated gene list files.\n'
						'Each file in this folder is named after a GO term (GO_#######) and each line of each file has a single gene name (Cc##_g#####)')
	parser.add_argument('output_folder', help='output folder for .gmt file and enrichments files')
	parser = parser.parse_args()

	# Handle folders
	if (os.path.exists(parser.genes_folder_up) and
		os.path.exists(parser.genes_folder_down) and
		os.path.exists(parser.output_folder)):
		
		# enrichments
		print('Parsing enrichments files ... ', end='')
		enrichments_header, enrichments_pos = enrichments(parser.enrichments_pos)  # .gmt header
		enrichments_header, enrichments_neg = enrichments(parser.enrichments_neg)  # {'GO:#######' : <gsea output line>, ... }
		print('OK')

		# go_and_description
		print('Parsing go_and_description files ... ', end='')
		go_and_description_up = go_and_description(parser.go_and_description_up)  # {'GO:#######' : 'GO:####### \t description', ... }
		go_and_description_down = go_and_description(parser.go_and_description_down)
		print('OK')

		# genes_folder
		print('Parsing genes folders ... ', end='')
		genes_folder_up = genes_folder(parser.genes_folder_up)  # {'GO:#######' : set('Cc##_g#####', ...), ... }
		genes_folder_down = genes_folder(parser.genes_folder_down)
		print('OK')

		# build gmt lines
		print('Building .gmt file lines ... ', end='')
		gmt_lines = {}  # {'GO:#######' : <.gmt line>, ... }
		go_and_description_merged = {**go_and_description_up, **go_and_description_down}
		for GO in go_and_description_merged:
			if GO in genes_folder_up:
				gmt_lines[GO] = go_and_description_merged[GO] + '\t' + '\t'.join(genes_folder_up[GO])
			if GO in genes_folder_down:
				if GO not in gmt_lines:
					gmt_lines[GO] = go_and_description_merged[GO] + '\t' + '\t'.join(genes_folder_down[GO])
				else:
					gmt_lines[GO] += '\t' + '\t'.join(genes_folder_up[GO])
		print('OK')
		
		# filter enrichments files for go terms
		print('Filter enrichments files for go terms that don\'t exist in the .gmt file ... ', end='')
		keys_to_delete_pos = []
		for GO in enrichments_pos:
			if GO not in gmt_lines:
				keys_to_delete_pos.append(GO)
		for key in keys_to_delete_pos:
			del enrichments_pos[key]
		
		keys_to_delete_neg = []
		for GO in enrichments_neg:
			if GO not in gmt_lines:
				keys_to_delete_neg.append(GO)
		for key in keys_to_delete_neg:
			del enrichments_neg[key]
		print('OK')

		# write output files
		print('Writing output files ... ', end='')
		with open(os.path.join(parser.output_folder, 'enrichmentMap.gmt'), 'w', newline='') as gmt:
			for gmt_line in gmt_lines.values():
				gmt.write(gmt_line + '\n')
		
		with open(os.path.join(parser.output_folder, os.path.basename(parser.enrichments_pos.name).split('.')[0] + '_enrichments_pos.txt'), 'w', newline='') as pos, \
			open(os.path.join(parser.output_folder, os.path.basename(parser.enrichments_neg.name).split('.')[0] + '_enrichments_neg.txt'), 'w', newline='') as neg:
			pos.write(enrichments_header)
			for line in enrichments_pos.values():
				pos.write(line + '\n')
			
			neg.write(enrichments_header)
			for line in enrichments_neg.values():
				neg.write(line + '\n')
		print('OK')

	else:
		print('Please input valid and existing paths for \'genes_folder_up\', \'genes_folder_down\' \'and output_folder\'.')
