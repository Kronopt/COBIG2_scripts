#!/usr/bin/env python3
# coding: utf-8

"""
Generates a table of the DEG's mean LogFC, up-regulated LogFC and down-regulated LogFC (filtered by significance, FDR > 0.1).
"""

import argparse
import csv
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generates a table of the DEG\'s mean LogFC, up-regulated LogFC and down-regulated LogFC (filtered by significance, FDR < 0.1).')
    parser.add_argument('de_results_folder', help='folder containing results.csv DE output files.')
    parser.add_argument('output_folder', help='output folder for .csv files.')
    parser = parser.parse_args()

    # Handle folders
    if os.path.exists(parser.de_results_folder) and os.path.exists(parser.output_folder):

        # de_results
        print('Parsing DE results.csv files ... ', end='')
        list_de_results = [de_results_file for de_results_file in os.listdir(parser.de_results_folder)
                           if 'results.csv' in de_results_file]
        de_results_dict = dict()
        for de_results_file in list_de_results:
            de_results_name = de_results_file.split('_')[0]
            de_results_dict[de_results_name] = {'up-sum': 0, 'up-total': 0, 'down-sum': 0, 'down-total': 0}
            with open(os.path.join(parser.de_results_folder, de_results_file), newline='') as de_file:
                de_file_reader = csv.reader(de_file)
                next(de_file_reader)  # ignore header
                for row in de_file_reader:
                    if row[6] != 'NA' and float(row[6]) < 0.1:  # FDR p-adjusted
                        if float(row[2]) < 0:  # log2FoldChange < 0
                            de_results_dict[de_results_name]['down-sum'] += float(row[2])
                            de_results_dict[de_results_name]['down-total'] += 1
                        elif float(row[2]) > 0:  # log2FoldChange > 0
                            de_results_dict[de_results_name]['up-sum'] += float(row[2])
                            de_results_dict[de_results_name]['up-total'] += 1
        print('OK')

        # DEG_mean_logfc.csv
        print('Generating DEG_mean_logfc.csv ... ', end='')
        with open(os.path.join(parser.output_folder, 'DEG_mean_logfc.csv'), 'w', newline='') as mean_logfc:
            values_csv_writer = csv.writer(mean_logfc)

            # header
            values_csv_writer.writerow(['Comparison', 'mean logFC', 'mean logFC up-regulated', 'mean logFC down-regulated'])

            for de_results_name in de_results_dict:  # for each comparison
                mean_up = de_results_dict[de_results_name]['up-sum'] / de_results_dict[de_results_name]['up-total']
                mean_down = de_results_dict[de_results_name]['down-sum'] / de_results_dict[de_results_name]['down-total']
                mean_total = ((de_results_dict[de_results_name]['up-sum'] + de_results_dict[de_results_name]['down-sum']) /
                             (de_results_dict[de_results_name]['up-total'] + de_results_dict[de_results_name]['down-total']))

                values_csv_writer.writerow([de_results_name, mean_total, mean_up, mean_down])
        print('OK')

    else:
        print('Please input valid and existing paths for \'de_results_folder\' and \'output_folder\'.')
