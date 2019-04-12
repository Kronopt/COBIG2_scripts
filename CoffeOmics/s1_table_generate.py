#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates S1 table (blast2go top10 go-terms).')
    parser.add_argument('seqtable_folder', help='folder containing seqtable.tab blast2go output files.')
    parser.add_argument('output_folder', help='output folder for S1 table.')
    parser = parser.parse_args()

    # Handle folders
    if (os.path.exists(parser.seqtable_folder) and
            os.path.exists(parser.output_folder)):

        # seqtables
        print('Parsing seqtable.tab files ... ', end='')
        list_seqtables = [seqtable_file for seqtable_file in os.listdir(parser.seqtable_folder)
                          if 'seqtable.tab' in seqtable_file]
        seqtable_go_count = dict()
        seqtable_go_count_total = dict()
        seqtable_go_count_top10 = dict()
        seqtable_total_gene_counts = set()

        for seqtable_file in list_seqtables:
            seqtable_name = seqtable_file.split('_')[0]  # ex: 1vs3
            up_or_down = 'down' if 'DOWN' in seqtable_file else 'up'
            if seqtable_name not in seqtable_go_count:
                seqtable_go_count[seqtable_name] = {'up': set(), 'down': set()}
            if seqtable_name not in seqtable_go_count_total:
                seqtable_go_count_total[seqtable_name] = {'C': dict(), 'F': dict(), 'P': dict()}
            if seqtable_name not in seqtable_go_count_top10:
                seqtable_go_count_top10[seqtable_name] = {'C': [], 'F': [], 'P': []}

            with open(os.path.join(parser.seqtable_folder, seqtable_file)) as seqtab_file:
                next(seqtab_file)  # ignore header
                for row in seqtab_file:
                    row = row.split('\t')
                    gene_name = row[0]
                    if len(row[6]) > 0 and int(row[6]) > 0:  # number of GO's > 0
                        seqtable_total_gene_counts.add(gene_name)
                        seqtable_go_count[seqtable_name][up_or_down].add(gene_name)

                        go_ids = row[7].split('; ')  # row[7] is GO ids
                        for go in go_ids:
                            go_domain = go[0]  # can be C, F or P
                            if go in seqtable_go_count_total[seqtable_name][go_domain]:
                                seqtable_go_count_total[seqtable_name][go_domain][go].add(gene_name)
                            else:
                                seqtable_go_count_total[seqtable_name][go_domain][go] = set()

            # top10, list of (key,value) pairs
            for go_domain in seqtable_go_count_total[seqtable_name]:
                seqtable_go_count_top10[seqtable_name][go_domain] = sorted(seqtable_go_count_total[seqtable_name][go_domain].items(), key=lambda z: len(z[1]), reverse=True)[:10]
        print('OK')

        print('Writing csv ... ', end='')
        # csv header
        header = ['GO-terms (top 10)', '# genes', '% total genes anotados', '% genes top 10',
                  '# genes down', '% total genes anotados', '% genes top 10',
                  '# genes up', '% total genes anotados', '% genes top 10']

        go_domain_name = {'C': 'Cell Component', 'F': 'Molecular Function', 'P': 'Biological Process'}
        for seqtable_name in seqtable_go_count_top10:
            with open(os.path.join(parser.output_folder, seqtable_name + '_S1_table.csv'), 'w', newline= '') as csv_file:
                csv_writer = csv.writer(csv_file)
                csv_writer.writerow(header)  # write header

                for go_domain in seqtable_go_count_top10[seqtable_name]:
                    csv_writer.writerow([go_domain_name[go_domain]])  # write GO top domain name

                    for go_term, gene_set in seqtable_go_count_top10[seqtable_name][go_domain]:
                        total_top10_genes = set()
                        for _, other_set in seqtable_go_count_top10[seqtable_name][go_domain]:
                            total_top10_genes = total_top10_genes.union(other_set)

                        genes_down = [gene_name for gene_name in gene_set
                                      if gene_name in seqtable_go_count[seqtable_name]['down']]
                        genes_up = [gene_name for gene_name in gene_set
                                      if gene_name in seqtable_go_count[seqtable_name]['up']]

                        row = []
                        row.append(go_term[2:])     # GO-terms (top 10)
                        row.append(str(         # # genes
                                   len(gene_set)))
                        row.append(str(         # % total genes anotados
                                   (len(gene_set) / len(seqtable_total_gene_counts)) * 100))
                        row.append(str(         # % genes top 10
                                   (len(gene_set) / len(total_top10_genes)) * 100))
                        row.append(str(         # # genes down
					               len(genes_down)))
                        row.append(str(         # % total genes anotados
                                   (len(genes_down) / len(seqtable_total_gene_counts)) * 100))
                        row.append(str(         # % genes top 10
                                   (len(genes_down) / len(total_top10_genes)) * 100))
                        row.append(str(         # # genes up
					               len(genes_up)))
                        row.append(str(         # % total genes anotados
                                   (len(genes_up) / len(seqtable_total_gene_counts)) * 100))
                        row.append(str(         # % genes top 10
                                   (len(genes_up) / len(total_top10_genes)) * 100))

                        csv_writer.writerow(row)

        print('OK')

    else:
        print('Please input valid and existing paths for \'seqtable_folder\' and \'output_folder\'.')
