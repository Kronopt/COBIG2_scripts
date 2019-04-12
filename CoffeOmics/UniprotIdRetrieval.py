#!python3
# coding: utf-8

"""
UNIPROT ID RETRIEVAL (adapted from UniprotIdRetrieval
Downloads protein sequences from Uniprot, in the desired output format, based on the given ids (either passed as
arguments or identified in a file)

UniprotIdRetrieval.py [-h] <file>

- '-h, --help' shows the help text
- 'file' takes a file path to search for ids
- '-fo, --format' defines the output format (defaults to fasta)
"""

import argparse
import re
import requests
import time


UNIPROTREGEX = ("[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]"
                + "|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9]")


def extract_uniprot_ids(file_name):
    uniprot_regex = re.compile(UNIPROTREGEX)
    uniprot_ids = []

    try:
        with open(file_name, "r") as file_with_ids:
            for line in file_with_ids:
                uniprot_ids.extend(uniprot_regex.findall(line))

            # Uniprot ids can only be 6 or 10 characters in size
            uniprot_ids = filter(lambda x: len(x) in [6, 10], uniprot_ids)

    except IOError as e:
        print('Error with "' + e.filename + '" file: ' + e.strerror)

    finally:
        return set(uniprot_ids)


def retrieve_sequences(ids_, output_format):
    if len(ids_) > 0:
        print('\n' + output_format.upper() + ' sequences')

        with open('log.txt', 'w') as log:
            for i in sorted(ids_):
                info_start = i + ': '
                print(info_start, end='')
                log.write(info_start)

                resource = i + '.' + output_format

                # Uniprot webservice
                sequence_file = requests.get('https://www.uniprot.org/uniprot/' + resource)

                # If response is empty
                if len(sequence_file.text) == 0:
                    seq_empty = 'not available in .' + output_format + ' or does not exist'
                    print(seq_empty)
                    log.write(seq_empty + '\n')
                    continue

                # http not 200
                elif sequence_file.status_code != 200:
                    seq_http_error = 'http error ' + str(sequence_file.status_code)
                    print(seq_http_error)
                    log.write(seq_http_error + '\n')

                # If response is html, then it's invalid
                else:
                    html = False
                    seq_html = 'not available in .' + output_format + ' or does not exist'
                    for line in sequence_file.iter_lines():
                        line = line.decode('utf-8')
                        if '<!DOCTYPE html' in line:
                            print(seq_html)
                            log.write(seq_html + '\n')
                            html = True
                        break

                    if html:
                        continue

                with open(resource, "w") as file_name:
                    [file_name.write(line.decode('utf-8') + '\n') for line in sequence_file.iter_lines()]

                print('ok')
                log.write('ok\n')
                time.sleep(0.33)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Uniprot id retrieval tool. Searches Uniprot for protein sequences and'
                                                 ' downloads them. Uses ids found in a file')

    parser.add_argument('f', metavar='<file>', help='file path')
    parser.add_argument('-fo', '--format',
                        choices=['fasta', 'gff', 'tab', 'txt', 'rdf', 'xml'],
                        default='fasta',
                        help='output format (defaults to fasta)')

    parser = parser.parse_args()

    ids = extract_uniprot_ids(parser.f)
    retrieve_sequences(ids, parser.format)
