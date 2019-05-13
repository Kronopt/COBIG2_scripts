#!/usr/bin/env python3
# coding: utf-8

"""
Use QuickGO API (https://www.ebi.ac.uk/QuickGO/) to get the descendants of GO terms.
Verify which genes are annotated with go terms descendant from those.
"""

import argparse, csv, json, os, requests, sys, urllib

# Photosynthesis				P:GO:0015979
# Chlorophyll					P:GO:0015994
# Pyruvate kinase				F:GO:0004743
# Malate dehydrogenase			F:GO:0016615
# Rubisco/Ribulose bisphosphate	F:GO:0016984
# get_descendants(['GO:0015979', 'GO:0015994', 'GO:0004743', 'GO:0016615', 'GO:0016984'])

URL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/descendants?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"

def get_descendants(go_term):
	"""
	go_term is a list/tuple of GO terms in the form of 'GO:#######'
	returns a dictionary where key=go_term, value=(go_term_description, [descendants, ...])
	"""
	requestURL = URL.format(urllib.parse.quote(','.join(go_term)))
	r = requests.get(requestURL, headers={"Accept" : "application/json"})

	if not r.ok:
		r.raise_for_status()
		sys.exit()

	final_dict = {}
	responseBody = json.loads(r.text)
	for result in responseBody['results']:
		final_dict[result['id']] = (result['name'], result['descendants'])
	return final_dict

def get_seqtable_go_terms(seqtable_csv):
	"""
	seqtable_csv is a csv.reader object
	returns a dictionary where key=go_term, value=[(gene_id, gene_description), ...]
	"""
	final_dict = {}

	next(seqtable_csv)  # ignore header
	for row in seqtable_csv:
		row = row.split('\t')
		if len(row[6]) > 0 and int(row[6]) > 0:  # number of GO's > 0
			gene_id_and_description = (row[0], row[1])

			go_ids = row[7].split('; ')  # row[7] is GO ids
			for go in go_ids:
				go = go[2:]  # remove starting C, F or P

				if go not in final_dict:
					final_dict[go] = [gene_id_and_description]
				else:
					final_dict[go].append(gene_id_and_description)
	return final_dict


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Generates a table of genes annotated with GO terms descendants of the given GO terms.')
	parser.add_argument('go_terms', type=lambda arg: arg.split(','),
		help='Comma separated GO terms (no spaces).')
	parser.add_argument('seqtable_folder', help='folder containing seqtable.tab blast2go output files.')
	parser.add_argument('de_results_folder', help='folder containing results.csv DE output files.')
	parser.add_argument('output_folder', help='output folder for table (.csv file)')
	parser = parser.parse_args()

	# Handle folders
	if (os.path.exists(parser.seqtable_folder) and
		os.path.exists(parser.de_results_folder) and
		os.path.exists(parser.output_folder)):

		# quickGO API call
		print('Getting descendants ... ', end='')
		descendants = get_descendants(parser.go_terms)
		print('OK')

		# seqtables
		print('Parsing seqtable.tab files ... ', end='')
		seqtab_go_terms = []
		list_seqtables = [seqtable_file for seqtable_file in os.listdir(parser.seqtable_folder)
							if 'seqtable.tab' in seqtable_file]
		for seqtable_file in list_seqtables:
			with open(os.path.join(parser.seqtable_folder, seqtable_file)) as seqtab_file:
				seqtab_go_terms.append(get_seqtable_go_terms(seqtab_file))
		print('OK')

		# de_results
		print('Parsing DE results.csv files ... ', end='')
		list_de_results = [de_results_file for de_results_file in os.listdir(parser.de_results_folder)
							if 'results.csv' in de_results_file]
		de_results_dict = dict()
		for de_results_file in list_de_results:
			de_results_name = de_results_file.split('_')[0]
			de_results_dict[de_results_name] = dict()
			with open(os.path.join(parser.de_results_folder, de_results_file), newline='') as de_file:
				de_file_reader = csv.reader(de_file)
				next(de_file_reader)  # ignore header
				for row in de_file_reader:
					if row[6] != 'NA' and float(row[6]) < 0.05:  # FDR p-adjusted
						de_results_dict[de_results_name][row[0]] = row[2]  # log2FoldChange
		print('OK')

		# compare seqtable go terms vs quickGO descendants
		print('Checking which go terms are descendants ... ', end='')
		verified_descendants = {}  # key=(go_term, go_term_description), value=set((descendant_id, descendant_description), ...)

		for seqtable_dict in seqtab_go_terms:  # dict
			for go_term_seqtable in seqtable_dict:  # key=go_term, value=[(gene_id, gene_description), ...]
				for go_term_quickgo in descendants:    # key=go_term, value=(go_term_description, [descendants, ...])
					if go_term_seqtable in descendants[go_term_quickgo][1]:
						if (go_term_quickgo, descendants[go_term_quickgo][0]) not in verified_descendants:
							verified_descendants[(go_term_quickgo, descendants[go_term_quickgo][0])] = set(seqtable_dict[go_term_seqtable])
						else:
							verified_descendants[(go_term_quickgo, descendants[go_term_quickgo][0])].update(seqtable_dict[go_term_seqtable])
		print('OK')

		# write output file
		print('Writing output file ... ', end='')
		header = ['GO term', 'GO term description', 'Gene ID (descendant)', 'Gene Description', 'Foldchange 1vs3','Foldchange 5vs7', 'Foldchange 1vs5', 'Foldchange 3vs7']

		with open(os.path.join(parser.output_folder, 'table5.csv'), 'w', newline= '') as csv_file:
			csv_writer = csv.writer(csv_file)
			csv_writer.writerow(header)  # write header

			for go_term_and_description in verified_descendants:
				for descend_id, descend_descrip in verified_descendants[go_term_and_description]:

					csv_writer.writerow([*(go_term_and_description), descend_id, descend_descrip,
					de_results_dict['1vs3'].get(descend_id, ''), de_results_dict['5vs7'].get(descend_id, ''),
					de_results_dict['1vs5'].get(descend_id, ''), de_results_dict['3vs7'].get(descend_id, '')])
		print('OK')

	else:
		print('Please input valid and existing paths for \'seqtable_folder\' and \'output_folder\'.')
