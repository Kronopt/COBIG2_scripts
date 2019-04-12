#!/usr/bin/env python3

"""
Builds 2 FASTA files (up-regulated and down-regulated genes)
with sequences that correspond to genes present in the output
of DESEQ2 analysis

results.csv (DESEQ2 output) is first filtered for padj<0.05,
then fasta sequences of those genes are extracted from
c_canephora_cds.fna

Inputs:
    - Folder containing <x>vs<y>_results.csv files (DESEQ2 analysis output)
    - Genome (c_canephora_cds.fna)

Output:
    - Fasta files with sequences of the genes with padj<0.05,
      each file containing eiter up- or down- regulated genes

call the script from the command line with '-h' to see documentation
"""


import argparse
import csv
import os


################
# HELPER CLASSES (adapted from PyFastaParser)
################


class FastaSequence:
    def __init__(self, seq_id, description, sequence):
        self.id = seq_id
        self.description = description
        self.sequence = sequence


class FastaParser:
    def __init__(self, fasta_file, keep_sequences=False):
        if isinstance(fasta_file, str):  # try to open file path
            self._fasta = open(fasta_file, 'r')
        elif hasattr(fasta_file, "readline"):  # assume it's a file object
            if fasta_file.closed:
                self._fasta = open(fasta_file.name, 'r')
            else:
                self._fasta = fasta_file
        else:
            raise Exception("Not a file.")

        self._iteration_ended = False
        self._keep_sequences = keep_sequences
        self.sequences = {}

    def __iter__(self):
        # check if file was closed, and open it again for new iteration
        if self._fasta.closed:
            self._fasta = open(self._fasta.name, 'r')

        def iter_fasta_file(fasta_file):
            # assumes first line begins with '>'
            first_line = fasta_file.readline()[1:].split(" ", 1)

            end_of_file = False
            while not end_of_file:
                if len(first_line) == 1:  # description can be empty
                    seq_id = first_line[0].rstrip()
                    seq_description = ''
                else:
                    seq_id, seq_description = first_line

                seq = ''
                end_of_sequence = False
                while not end_of_sequence:
                    sequence_line = fasta_file.readline()
                    if len(sequence_line) == 0:  # end of file, end iteration
                        end_of_sequence = True
                        end_of_file = True
                        fasta_file.close()
                    elif not sequence_line.startswith('>'):  # Line containing part of the sequence
                        seq += sequence_line.rstrip()
                    else:
                        end_of_sequence = True
                        first_line = sequence_line[1:].split(" ", 1)

                fasta_sequence = FastaSequence(seq_id, seq_description.rstrip(), seq)

                if self._keep_sequences and not self._iteration_ended:
                    self.sequences[fasta_sequence.id] = fasta_sequence

                yield fasta_sequence

            self._iteration_ended = True  # Iterated once. No more sequences will be added to self.sequences

        return iter_fasta_file(self._fasta)


def main(input_csv, input_genome, output_file_up, output_file_down):
    input_csv_reader = csv.reader(input_csv, quoting=csv.QUOTE_NONNUMERIC)

    # ignore header on csv input
    next(input_csv_reader)

    # read and filter csv input in memory
    # dictionary for faster input_csv_filtered_rows appending
    input_csv_filtered_rows = []
    gene_name_position = {}
    count = 0
    for line in input_csv_reader:
        if float(line[6]) < 0.05:
            input_csv_filtered_rows.append(line)
            gene_name_position[line[0]] = count
            count += 1

    # for each gene_name, add that gene fasta sequence (with header) to the output file
    for line in input_csv_filtered_rows:
        gene_name = line[0]
        fasta_sequence = input_genome.sequences[gene_name]

        # up- vs down- regulated (log2FoldChange)
        if line[2] > 0:
            output_file_up.write('>' + fasta_sequence.id + ' ' + fasta_sequence.description + '\n')
            output_file_up.write(fasta_sequence.sequence + '\n')
        elif line[2] < 0:
            output_file_down.write('>' + fasta_sequence.id + ' ' + fasta_sequence.description + '\n')
            output_file_down.write(fasta_sequence.sequence + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Builds 2 FASTA files (up- and down- regulated genes) with sequences that correspond to genes present in the output of DESEQ2 analysis.')
    parser.add_argument('results_csv_folder', help='Folder containing the input results.csv files (DESEQ2 output).')
    parser.add_argument('genome', type=argparse.FileType('r', encoding='UTF-8'),
                        help='genome file (c_canephora_cds.fna).')
    parser.add_argument('output_folder', help='output folder for fasta files.')
    parser = parser.parse_args()

    # Parse whole genome to memory
    print('Parsing whole genome file to memory ... ', end='')
    fasta_parser = FastaParser(parser.genome, keep_sequences=True)
    for seq in fasta_parser:
        pass
    print('OK')

    # If folders exist
    if os.path.exists(parser.results_csv_folder) and os.path.exists(parser.output_folder):
        list_of_results_csv = [csv_file for csv_file in os.listdir(parser.results_csv_folder) if csv_file.endswith('results.csv')]

        for results_csv in list_of_results_csv:
            results_csv_path_input = os.path.join(parser.results_csv_folder, results_csv)
            results_csv_path_output_up = os.path.join(parser.output_folder, ''.join([results_csv[:results_csv.find('.csv')], '_UP.fasta']))
            results_csv_path_output_down = os.path.join(parser.output_folder, ''.join([results_csv[:results_csv.find('.csv')], '_DOWN.fasta']))

            with open(results_csv_path_input, encoding='UTF-8') as csv_input, open(results_csv_path_output_up, 'w', encoding='UTF-8') as fasta_output_up, open(results_csv_path_output_down, 'w', encoding='UTF-8') as fasta_output_down:
                print('parsing', results_csv, '... ', end='')
                main(csv_input, fasta_parser, fasta_output_up, fasta_output_down)
                print('OK')

    else:
        print('Please input valid and existing paths for \'results_csv_folder\' and \'output_folder\'')
