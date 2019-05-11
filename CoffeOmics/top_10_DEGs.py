#!/usr/bin/env python3
# coding: utf-8

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
        for de_results_name in de_results_dict:
            top_10_de_up_list = sorted(de_results_dict[de_results_name]['up'].items(), key=lambda z: z[1], reverse=True)[:10]
            top_10_de_down_list = sorted(de_results_dict[de_results_name]['down'].items(), key=lambda z: z[1], reverse=True)[:10]

            with open(os.path.join(parser.output_folder, de_results_name + '_top10_DEGs.csv'), 'w', newline='') as top10_csv:
                values_csv_writer = csv.writer(top10_csv)

                # up-regulated genes
                values_csv_writer.writerow([de_results_name + ' UP-regulated'])  # sub-header 
                for gene_name, _ in top_10_de_up_list:
                    # genes_description = seqtable_dict[de_results_name]['up'][gene_name]
                    try:
                        genes_description = seqtable_dict[gene_name]
                    except KeyError:
                        genes_description = ''
                    values_csv_writer.writerow([gene_name, genes_description])

                # down-regulated genes
                values_csv_writer.writerow([de_results_name + ' DOWN-regulated'])  # sub-header
                for gene_name, _ in top_10_de_down_list:
                    try:
                        genes_description = seqtable_dict[gene_name]
                    except KeyError:
                        genes_description = ''
                    values_csv_writer.writerow([gene_name, genes_description])
        print('OK')

    else:
        print('Please input valid and existing paths for \'de_results_folder\',',
              '\'seqtable_folder\' and \'output_folder\'.')
