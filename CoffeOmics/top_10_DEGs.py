#!/usr/bin/env python3
# coding: utf-8

"""
Generates a table of the top 10 DEG's UP- and DOWN-regulated for each comparison.
Also generates the top 10 DEG's UP- and DOWN-regulated for 1v5 VS 3v7 and 1v3 VS 5v7 specific to each comparison
(top 10 DEG's of 1v5 not in 3v7 and vice versa, same for 1v3 VS 5v7)
"""

import argparse
import csv
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generates a table of top 10 DEGs (differential expressed genes).')
    parser.add_argument('de_results_folder', help='folder containing results.csv DE output files.')
    parser.add_argument('seqtable_folder', help='folder containing seqtable.tab blast2go output files.')
    parser.add_argument('output_folder', help='output folder for .csv files.')
    parser = parser.parse_args()

    # Handle folders
    if (os.path.exists(parser.de_results_folder) and
            os.path.exists(parser.seqtable_folder) and
            os.path.exists(parser.output_folder)):

        # de_results
        print('Parsing DE results.csv files ... ', end='')
        list_de_results = [de_results_file for de_results_file in os.listdir(parser.de_results_folder)
                           if 'results.csv' in de_results_file]
        de_results_dict = dict()
        for de_results_file in list_de_results:
            de_results_name = de_results_file.split('_')[0]
            de_results_dict[de_results_name] = {'up': dict(), 'down': dict()}
            with open(os.path.join(parser.de_results_folder, de_results_file), newline='') as de_file:
                de_file_reader = csv.reader(de_file)
                next(de_file_reader)  # ignore header
                for row in de_file_reader:
                    if row[6] != 'NA' and float(row[6]) < 0.05:  # FDR p-adjusted
                        if float(row[2]) < 0:  # log2FoldChange < 0
                            de_results_dict[de_results_name]['down'][row[0]] = row[2]
                        elif float(row[2]) > 0:  # log2FoldChange > 0
                            de_results_dict[de_results_name]['up'][row[0]] = row[2]
        print('OK')

        # seqtables
        print('Parsing seqtable.tab files ... ', end='')
        list_seqtables = [seqtable_file for seqtable_file in os.listdir(parser.seqtable_folder)
                          if 'seqtable.tab' in seqtable_file]
        seqtable_dict = dict()
        for seqtable_file in list_seqtables:
            # seqtable_name = seqtable_file.split('_')[0]
            # up_or_down = 'down' if 'DOWN' in seqtable_file else 'up'
            # if seqtable_name not in seqtable_dict:
            #     seqtable_dict[seqtable_name] = {'up': dict(), 'down': dict()}
            with open(os.path.join(parser.seqtable_folder, seqtable_file)) as seqtab_file:
                next(seqtab_file)  # ignore header
                for row in seqtab_file:
                    row = row.split('\t')
                    if row[1] != '---NA---':  # has a sequence description
                        # seqtable_dict[seqtable_name][up_or_down][row[0]] = row[1]
                        seqtable_dict[row[0]] = row[1]
        print('OK')

        # top10_DEGs.csv
        print('Generating top10_DEGs.csv ... ', end='')
        ordered_de = {}
        for de_results_name in de_results_dict:
            ordered_de[de_results_name] = {}
            ordered_de[de_results_name]['up'] = sorted(de_results_dict[de_results_name]['up'].items(), key=lambda z: z[1], reverse=True)
            ordered_de[de_results_name]['down'] = sorted(de_results_dict[de_results_name]['down'].items(), key=lambda z: z[1], reverse=True)

            top_10_de_up_list = ordered_de[de_results_name]['up'][:10]
            top_10_de_down_list = ordered_de[de_results_name]['down'][:10]

            with open(os.path.join(parser.output_folder, de_results_name + '_top10_DEGs.csv'), 'w', newline='') as top10_csv:
                values_csv_writer = csv.writer(top10_csv)

                # up-regulated genes
                values_csv_writer.writerow([de_results_name + ' UP-regulated'])  # sub-header 
                for gene_name, foldchange in top_10_de_up_list:
                    # genes_description = seqtable_dict[de_results_name]['up'][gene_name]
                    try:
                        genes_description = seqtable_dict[gene_name]
                    except KeyError:
                        genes_description = ''
                    values_csv_writer.writerow([gene_name, genes_description, foldchange])

                # down-regulated genes
                values_csv_writer.writerow([de_results_name + ' DOWN-regulated'])  # sub-header
                for gene_name, foldchange in top_10_de_down_list:
                    try:
                        genes_description = seqtable_dict[gene_name]
                    except KeyError:
                        genes_description = ''
                    values_csv_writer.writerow([gene_name, genes_description, foldchange])
        print('OK')

        # top10_DEGs_comparison_specific.csv
        print('Generating top10_DEGs_comparison_specific.csv ... ', end='')
        comparions_to_do = {'1vs5': '3vs7', '3vs7': '1vs5', '1vs3': '5vs7', '5vs7': '1vs3'}
        comparisons_names = {'1vs5': '1vs5 VS 3vs7', '3vs7': '1vs5 VS 3vs7', '1vs3': '1vs3 vs 5vs7', '5vs7': '1vs3 vs 5vs7'}
        comparisons = {'1vs5 VS 3vs7': {
                          'up': {'1vs5': [], '3vs7': []},
                          'down': {'1vs5': [], '3vs7': []}
                      },'1vs3 vs 5vs7': {
                          'up': {'1vs3': [], '5vs7': []},
                          'down': {'1vs3': [], '5vs7': []}}}
        
        # For each comparison (1vs5, 3vs7, 1vs3, 5vs7) check the top genes that are not in the other comparison of the same up/down regulated set.
        # Do this until 10 genes are found
        for de_results_name in de_results_dict:
            gene_count = 10
            for gene_name, foldchange in ordered_de[de_results_name]['up']:
                if gene_name not in de_results_dict[comparions_to_do[de_results_name]]['up']:
                    try:
                        genes_description = seqtable_dict[gene_name]
                    except KeyError:
                        genes_description = ''
                    comparisons[comparisons_names[de_results_name]]['up'][de_results_name].append((gene_name, genes_description, foldchange))
                    gene_count -= 1
                    if gene_count == 0:
                        break
            gene_count = 10
            for gene_name, foldchange in ordered_de[de_results_name]['down']:
                if gene_name not in de_results_dict[comparions_to_do[de_results_name]]['down']:
                    try:
                        genes_description = seqtable_dict[gene_name]
                    except KeyError:
                        genes_description = ''
                    comparisons[comparisons_names[de_results_name]]['down'][de_results_name].append((gene_name, genes_description, foldchange))
                    gene_count -= 1
                    if gene_count == 0:
                        break
        for comparison_of_comparisons in comparisons:
            with open(os.path.join(parser.output_folder, comparison_of_comparisons.replace(' ', '_') + '_top10_DEGs.csv'), 'w', newline='') as top10_csv:
                values_csv_writer = csv.writer(top10_csv)
                values_csv_writer.writerow(['UP-regulated'])  # sub-header
                for comp in comparisons[comparison_of_comparisons]['up']:
                    values_csv_writer.writerow([comp + ' only'])  # sub-header 2
                    values_csv_writer.writerows(comparisons[comparison_of_comparisons]['up'][comp])
                values_csv_writer.writerow(['DOWN-regulated'])  # sub-header
                for comp in comparisons[comparison_of_comparisons]['down']:
                    values_csv_writer.writerow([comp + ' only'])  # sub-header 2
                    values_csv_writer.writerows(comparisons[comparison_of_comparisons]['down'][comp])
        print('OK')

    else:
        print('Please input valid and existing paths for \'de_results_folder\',',
              '\'seqtable_folder\' and \'output_folder\'.')
