#!/usr/bin/env python3

"""
Extracts SwissProt IDs of each gene from a gff3 functional annotation file.

Inputs:
    - c_canephora_gene_structural_and_functional_annotation.gff3

Outputs:
    - .txt file listing SwissProt IDs
      (every SwissProt ID available for each gene)

call the script from the command line with '-h' to see documentation
"""


import argparse


def main(input_gff3, output_file):
    uniprotids = set()

    # fetch SwissProt IDs
    for line in input_gff3:
        # check if line has Swissprot IDs
        db_ref = line.find('Dbxref=')
        if db_ref != -1:
            
            swissprotids = []
            there_are_more_swissprot_ids = True
            line_location = 0
            while there_are_more_swissprot_ids:
                swissprot_id = line.find('SwissProt:', line_location)
                if swissprot_id != -1:
                    end_of_swissprot_id = line.find(',', swissprot_id)
                    if end_of_swissprot_id - swissprot_id > 20:  # the end of an id might be ';'
                        end_of_swissprot_id = line.find(';', swissprot_id)
                    swissprot_id_value = line[swissprot_id + 10: end_of_swissprot_id]
                    line_location = end_of_swissprot_id
                    swissprotids.append(swissprot_id_value)
                else:
                    there_are_more_swissprot_ids = False

            for swissprotid in swissprotids:
                uniprotids.add(swissprotid)

    # write everything to the output csv files
    for uniprotid in uniprotids:
        output_file.write(uniprotid + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extracts SwissProt IDs of each gene from a gff3 functional annotation file.')
    parser.add_argument('gff3', type=argparse.FileType('r', encoding='UTF-8'),
                        help='input gff3 file (genome functional annotation).')
    parser.add_argument('output_file', type=argparse.FileType('w', encoding='UTF-8'),
                        help='.txt file containing every Swissprot ID available for each gene in the gff3 file.')
    parser = parser.parse_args()

    main(parser.gff3, parser.output_file)
