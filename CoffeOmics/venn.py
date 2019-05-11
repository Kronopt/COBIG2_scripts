#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import os
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot


##############
# FASTA PARSER
##############


nucleotide_letter_codes_good = {
    'A': 'adenosine',
    'C': 'cytidine',
    'G': 'guanine',
    'T': 'thymidine',
    'N': 'any (A/G/C/T)',
    'U': 'uridine'
}
nucleotide_letter_codes_degenerate = {
    'K': 'keto (G/T)',
    'S': 'strong (G/C)',
    'Y': 'pyrimidine (T/C)',
    'M': 'amino (A/C)',
    'W': 'weak (A/T)',
    'R': 'purine (G/A)',
    'B': 'G/T/C',
    'D': 'G/A/T',
    'H': 'A/C/T',
    'V': 'G/C/A',
    '-': 'gap of indeterminate length'
}

aminoacid_letter_codes_good = {
    'A': 'alanine',
    'B': 'aspartate/asparagine',
    'C': 'cystine',
    'D': 'aspartate',
    'E': 'glutamate',
    'F': 'phenylalanine',
    'G': 'glycine',
    'H': 'histidine',
    'I': 'isoleucine',
    'K': 'lysine',
    'L': 'leucine',
    'M': 'methionine',
    'N': 'asparagine',
    'P': 'proline',
    'Q': 'glutamine',
    'R': 'arginine',
    'S': 'serine',
    'T': 'threonine',
    'U': 'selenocysteine',
    'V': 'valine',
    'W': 'tryptophan',
    'Y': 'tyrosine',
    'Z': 'glutamate/glutamine',
    'X': 'any',
    '*': 'translation stop'
}
aminoacid_letter_codes_degenerate = {
    '-': 'gap of indeterminate length'
}

# set operations
nucleotide_letter_codes_all = set(list(nucleotide_letter_codes_good) +
                                  list(nucleotide_letter_codes_degenerate))
aminoacid_letter_codes_all = set(list(aminoacid_letter_codes_good) +
                                 list(aminoacid_letter_codes_degenerate))
aminoacids_not_in_nucleotides = aminoacid_letter_codes_all - nucleotide_letter_codes_all
# nucleotides_not_in_aminoacids would be empty


class LetterCode:
    _dictionary = {
        'nucleotide': (nucleotide_letter_codes_good, nucleotide_letter_codes_degenerate),
        'aminoacid': (aminoacid_letter_codes_good, aminoacid_letter_codes_degenerate)
    }

    def __init__(self, letter_code, sequence_type=''):
        self.letter_code = letter_code.upper()
        self.description = ''
        self.degenerate = None
        self.supported = True

        if sequence_type:
            # letter_codes_good
            if self.letter_code in self._dictionary[sequence_type][0]:
                self.description = self._dictionary[sequence_type][0][self.letter_code]
                self.degenerate = False
            # letter_codes_degenerate
            elif self.letter_code in self._dictionary[sequence_type][1]:
                self.description = self._dictionary[sequence_type][1][self.letter_code]
                self.degenerate = True
            else:
                self.supported = False
        else:
            self.supported = False

    def __repr__(self):
        return self.letter_code


class FastaSequence:
    def __init__(self, seq_id, description, sequence, sequence_type='', infer_type=False):
        self.id = seq_id
        self.description = description
        self.sequence = sequence
        self.sequence_type = sequence_type
        self._infer_type = infer_type

        if self._infer_type:
            self._infer_sequence_type()

    def _infer_sequence_type(self):
        for letter_code in self.sequence:
            if letter_code in aminoacids_not_in_nucleotides:
                self.sequence_type = 'aminoacid'
                return
        self.sequence_type = ''

    def __iter__(self):

        def iter_sequence(sequence, sequence_type):
            for character in sequence:
                yield LetterCode(character, sequence_type)

        return iter_sequence(self.sequence, self.sequence_type)

    def __len__(self):
        return len(self.sequence)

    def __repr__(self):
        _id = self.id if self.id else '\'\''
        description = self.description[:50] if self.description else '\'\''
        return "<%s - ID:%s | DESCRIPTION:%s>" % (self.__class__.__name__, _id, description)


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

    def __repr__(self):
        return "<%s - FASTAFILE:%s>" % (self.__class__.__name__, self._fasta.name)


######
# MAIN
######


def create_venn(list_of_sets, list_of_labels, output_folder, _format):
    if len(list_of_sets) == 2:
        venn2(list_of_sets, list_of_labels)
        # values.csv
        intersect = list_of_sets[0].intersection(list_of_sets[1])
        venn_values.append([
            ' vs '.join(list_of_labels),
            str(len(list_of_sets[0]) - len(intersect)),
            str(len(list_of_sets[0])),
            str(len(intersect)),
            str(len(list_of_sets[1])),
            str(len(list_of_sets[1]) - len(intersect))
        ])
    elif len(list_of_sets) == 3:
        venn3(list_of_sets, list_of_labels)
    else:
        raise(ValueError, '\'sets\' can either be of length 2 or 3.')
    title = pyplot.title(' vs '.join(list_of_labels))
    for format_type in _format:
        pyplot.savefig(os.path.join(output_folder, title._text.replace(' ', '_').lower() + '.' + format_type), format=format_type)
    pyplot.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Creates venn diagrams for total gene expression, differential gene expression and GO annotated genes.')
    parser.add_argument('genome', type=argparse.FileType('r', encoding='UTF-8'),
                        help='genome file (c_canephora_cds.fna).')
    parser.add_argument('htseq_counts_folder', help='folder containing htseq-counts files.')
    parser.add_argument('de_results_folder', help='folder containing results.csv DE output files.')
    parser.add_argument('seqtable_folder', help='folder containing seqtable.tab blast2go output files.')
    parser.add_argument('output_folder', help='output folder for venn diagrams.')
    parser = parser.parse_args()

    venn_values = []  # for saving values used to generate venns
    venn_values.append(['Name', 'only A', 'total A', 'intersect','total B', 'only B'])  # header

    # Parse whole genome to memory
    print('Parsing whole genome file to memory ... ', end='')
    fasta_parser = FastaParser(parser.genome, keep_sequences=True)
    for seq in fasta_parser:
        pass
    genome_set = set(fasta_parser.sequences.keys())
    print('OK')

    # Handle folders
    if (os.path.exists(parser.htseq_counts_folder) and
            os.path.exists(parser.de_results_folder) and
            os.path.exists(parser.seqtable_folder) and
            os.path.exists(parser.output_folder)):

        # htseq_counts
        print('Parsing htseq-count files ... ', end='')
        list_htseq_counts = [htseq_count_file for htseq_count_file in os.listdir(parser.htseq_counts_folder)
                             if 'htseq_counts' in htseq_count_file]
        htseq_counts_set_dict = dict()
        for htseq_count_file in list_htseq_counts:
            # htseq-count file name: trimmed_<#><A>(...) or trimmed_<##><A>(...)
            if htseq_count_file[9] in ('A', 'B', 'C'):
                htseq_name = int(htseq_count_file[8:9])
            else:
                htseq_name = int(htseq_count_file[8:10])
            if htseq_name not in htseq_counts_set_dict:
                htseq_counts_set_dict[htseq_name] = set()
            with open(os.path.join(parser.htseq_counts_folder, htseq_count_file), newline='') as htseq_file:
                htseq_count_file_reader = csv.reader(htseq_file, delimiter='\t')
                for row in htseq_count_file_reader:
                    if int(row[1]) > 0:  # count > 0
                        htseq_counts_set_dict[htseq_name].add(row[0])
        print('OK')

        # de_results
        print('Parsing DE results.csv files ... ', end='')
        list_de_results = [de_results_file for de_results_file in os.listdir(parser.de_results_folder)
                           if 'results.csv' in de_results_file]
        de_results_set_dict = dict()
        for de_results_file in list_de_results:
            de_results_name = de_results_file[:-12]
            de_results_set_dict[de_results_name] = {'up': set(), 'down': set()}
            with open(os.path.join(parser.de_results_folder, de_results_file), newline='') as de_file:
                de_file_reader = csv.reader(de_file)
                next(de_file_reader)  # ignore header
                for row in de_file_reader:
                    if row[6] != 'NA' and float(row[6]) < 0.05:  # FDR p-adjusted
                        if float(row[2]) < 0:  # log2FoldChange < 0
                            de_results_set_dict[de_results_name]['down'].add(row[0])
                        elif float(row[2]) > 0:  # log2FoldChange > 0
                            de_results_set_dict[de_results_name]['up'].add(row[0])
        print('OK')

        # seqtables
        print('Parsing seqtable.tab files ... ', end='')
        list_seqtables = [seqtable_file for seqtable_file in os.listdir(parser.seqtable_folder)
                          if 'seqtable.tab' in seqtable_file]
        seqtable_set_dict = dict()
        for seqtable_file in list_seqtables:
            seqtable_name = seqtable_file.split('_')[0]
            up_or_down = 'down' if 'DOWN' in seqtable_file else 'up'
            if seqtable_name not in seqtable_set_dict:
                seqtable_set_dict[seqtable_name] = {'up': set(), 'down': set()}
            with open(os.path.join(parser.seqtable_folder, seqtable_file)) as seqtab_file:
                next(seqtab_file)  # ignore header
                for row in seqtab_file:
                    row = row.split('\t')
                    if len(row[6]) > 0 and int(row[6]) > 0:  # number of GO's
                        seqtable_set_dict[seqtable_name][up_or_down].add(row[0])
        print('OK')

        #######
        # VENNS
        #######

        print('Generating venn diagrams ... ', end='')

        # # c.canephora genome VS all expressed genes:
        # # all htseq-counts (count > 0)
        # create_venn([genome_set,
                     # set.union(*[htseq_set for htseq_set in htseq_counts_set_dict.values()])],
                    # ['c.canephora genome', 'all expressed genes'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS C380:
        # # htseq_counts 1, 2, 5, 6, 9, 11 (count > 0)
        # create_venn([genome_set,
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 2, 5, 6, 9, 11]])],
                    # ['c.canephora genome', 'C380'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS C700:
        # # htseq_counts 3, 4, 7, 8, 10, 12 (count > 0)
        # create_venn([genome_set,
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3, 4, 7, 8, 10, 12]])],
                    # ['c.canephora genome', 'C700'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS C380 VS C700:
        # # htseq_counts 1, 2, 5, 6, 9, 11 (count > 0)
        # # htseq_counts 3, 4, 7, 8, 10, 12 (count > 0)
        # create_venn([genome_set,
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 2, 5, 6, 9, 11]]),
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3, 4, 7, 8, 10, 12]])],
                    # ['c.canephora genome', 'C380', 'C700'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # C380 VS C700:
        # # htseq_counts 1, 2, 5, 6, 9, 11 (count > 0)
        # # htseq_counts 3, 4, 7, 8, 10, 12 (count > 0)
        # create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 2, 5, 6, 9, 11]]),
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3, 4, 7, 8, 10, 12]])],
                    # ['C380', 'C700'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS 25C:
        # # htseq_counts 1, 3, 5, 7 (count > 0)
        # create_venn([genome_set,
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 3, 5, 7]])],
                    # ['c.canephora genome', '25C'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS HiC:
        # # htseq_counts 2, 4, 6, 8 (count > 0)
        # create_venn([genome_set,
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [2, 4, 6, 8]])],
                    # ['c.canephora genome', 'HiC'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS ItC:
        # # htseq_counts 9, 10, 11, 12 (count > 0)
        # create_venn([genome_set,
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [9, 10, 11, 12]])],
                    # ['c.canephora genome', 'ItC'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # 25C VS HiC:
        # # htseq_counts 1, 3, 5, 7 (count > 0)
        # # htseq_counts 2, 4, 6, 8 (count > 0)
        # create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 3, 5, 7]]),
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [2, 4, 6, 8]])],
                    # ['25C', 'HiC'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # 25C VS ItC:
        # # htseq_counts 1, 3, 5, 7 (count > 0)
        # # htseq_counts 9, 10, 11, 12 (count > 0)
        # create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 3, 5, 7]]),
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [9, 10, 11, 12]])],
                    # ['25C', 'ItC'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # HiC VS ItC:
        # # htseq_counts 2, 4, 6, 8 (count > 0)
        # # htseq_counts 9, 10, 11, 12 (count > 0)
        # create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [2, 4, 6, 8]]),
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [9, 10, 11, 12]])],
                    # ['HiC', 'ItC'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # 25C VS HiC VS ItC:
        # # htseq_counts 1, 3, 5, 7 (count > 0)
        # # htseq_counts 2, 4, 6, 8 (count > 0)
        # # htseq_counts 9, 10, 11, 12 (count > 0)
        # create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 3, 5, 7]]),
                     # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [2, 4, 6, 8]]),
                    # set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [9, 10, 11, 12]])],
                    # ['25C', 'HiC', 'ItC'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora genome VS all DE genes:
        # # all results.csv (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([genome_set,
                     # set.union(*[de_results['up'].union(de_results['down'])
                                 # for de_results in de_results_set_dict.values()])],
                    # ['c.canephora genome', 'all DE genes'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # ######
        # # TEMP
        # ######

        # # c.canephora DE VS c.arabica DE:
        # # results.csv 5vs6, 5vs11, 6vs11, 7vs8, 7vs12, 8vs12 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 1vs2, 1vs9, 2vs9, 3vs4, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '6vs11', '7vs8', '7vs12', '8vs12']]),
                     # set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['1vs2', '1vs9', '2vs9', '3vs4', '3vs10', '4vs10']])],
                    # ['c.canephora DE', 'c.arabica DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora UP-DE VS c.arabica UP-DE:
        # # results.csv 5vs6, 5vs11, 6vs11, 7vs8, 7vs12, 8vs12 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 1vs2, 1vs9, 2vs9, 3vs4, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '6vs11', '7vs8', '7vs12', '8vs12']]),
                     # set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['1vs2', '1vs9', '2vs9', '3vs4', '3vs10', '4vs10']])],
                    # ['c.canephora up-DE', 'c.arabica up-DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # c.canephora DOWN-DE VS c.arabica DOWN-DE:
        # # results.csv 5vs6, 5vs11, 6vs11, 7vs8, 7vs12, 8vs12 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 1vs2, 1vs9, 2vs9, 3vs4, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '6vs11', '7vs8', '7vs12', '8vs12']]),
                     # set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['1vs2', '1vs9', '2vs9', '3vs4', '3vs10', '4vs10']])],
                    # ['c.canephora down-DE', 'c.arabica down-DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # C380 DE VS C700 DE:
        # # results.csv 5vs6, 5vs11, 6vs11, 1vs2, 1vs9, 2vs9 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 7vs8, 7vs12, 8vs12, 3vs4, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '6vs11', '1vs2', '1vs9', '2vs9']]),
                     # set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['7vs8', '7vs12', '8vs12', '3vs4', '3vs10', '4vs10']])],
                    # ['C380 DE', 'C700 DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # C380 UP-DE VS C700 UP-DE:
        # # results.csv 5vs6, 5vs11, 6vs11, 1vs2, 1vs9, 2vs9 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 7vs8, 7vs12, 8vs12, 3vs4, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '6vs11', '1vs2', '1vs9', '2vs9']]),
                     # set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['7vs8', '7vs12', '8vs12', '3vs4', '3vs10', '4vs10']])],
                    # ['C380 up-DE', 'C700 up-DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # C380 DOWN-DE VS C700 DOWN-DE:
        # # results.csv 5vs6, 5vs11, 6vs11, 1vs2, 1vs9, 2vs9 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 7vs8, 7vs12, 8vs12, 3vs4, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '6vs11', '1vs2', '1vs9', '2vs9']]),
                     # set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['7vs8', '7vs12', '8vs12', '3vs4', '3vs10', '4vs10']])],
                    # ['C380 down-DE', 'C700 down-DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # 25C DE VS HiC DE VS ItC DE:
        # # results.csv 5vs6, 5vs11, 1vs2, 1vs9, 7vs8, 7vs12, 3vs4, 3vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 5vs6, 6vs11, 1vs2, 2vs9, 7vs8, 8vs12, 3vs4, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 5vs11, 6vs11, 1vs9, 2vs9, 7vs12, 8vs12, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '1vs2', '1vs9', '7vs8', '7vs12', '3vs4', '3vs10']]),
                     # set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '6vs11', '1vs2', '2vs9', '7vs8', '8vs12', '3vs4', '4vs10']]),
                     # set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs11', '6vs11', '1vs9', '2vs9', '7vs12', '8vs12', '3vs10', '4vs10']])],
                    # ['25C DE', 'HiC DE', 'ItC DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # 25C UP-DE VS HiC UP-DE VS ItC UP-DE:
        # # results.csv 5vs6, 5vs11, 1vs2, 1vs9, 7vs8, 7vs12, 3vs4, 3vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 5vs6, 6vs11, 1vs2, 2vs9, 7vs8, 8vs12, 3vs4, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 5vs11, 6vs11, 1vs9, 2vs9, 7vs12, 8vs12, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '1vs2', '1vs9', '7vs8', '7vs12', '3vs4', '3vs10']]),
                     # set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '6vs11', '1vs2', '2vs9', '7vs8', '8vs12', '3vs4', '4vs10']]),
                     # set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs11', '6vs11', '1vs9', '2vs9', '7vs12', '8vs12', '3vs10', '4vs10']])],
                    # ['25C up-DE', 'HiC up-DE', 'ItC up-DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        # # 25C DOWN-DE VS HiC DOWN-DE VS ItC DOWN-DE:
        # # results.csv 5vs6, 5vs11, 1vs2, 1vs9, 7vs8, 7vs12, 3vs4, 3vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 5vs6, 6vs11, 1vs2, 2vs9, 7vs8, 8vs12, 3vs4, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # # results.csv 5vs11, 6vs11, 1vs9, 2vs9, 7vs12, 8vs12, 3vs10, 4vs10 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '5vs11', '1vs2', '1vs9', '7vs8', '7vs12', '3vs4', '3vs10']]),
                     # set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs6', '6vs11', '1vs2', '2vs9', '7vs8', '8vs12', '3vs4', '4vs10']]),
                     # set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 # if key in ['5vs11', '6vs11', '1vs9', '2vs9', '7vs12', '8vs12', '3vs10', '4vs10']])],
                    # ['25C down-DE', 'HiC down-DE', 'ItC down-DE'],
                    # parser.output_folder,
                    # ['png', 'eps'])

        #####
        # CO2
        #####

        # c.canephora VS c.arabica:
        # htseq_counts 5, 7 (count > 0)
        # htseq_counts 1, 3 (count > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [5, 7]]),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 3]])],
                    ['c.canephora expressed genes (5and7)', 'c.arabica expressed genes (1and3)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora GO VS c.arabica GO:
        # htseq_counts 5, 7 (count > 0)
        # htseq_counts 1, 3 (count > 0)
        # GO ALL (GO > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [5, 7]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()])),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 3]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()]))],
                    ['c.canephora expressed genes GO (5and7)', 'c.arabica expressed genes GO (1and3)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora DE VS c.arabica DE (co2 effect):
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]),
                     set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']])],
                    ['c.canephora DE (5vs7)(co2 effect)', 'c.arabica DE (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora DE GO VS c.arabica DE GO (co2 effect):
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 5vs7, 1vs3 (GO > 0)
        create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items() if key in ['5vs7']])),
                     set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items() if key in ['1vs3']]))],
                    ['c.canephora DE GO (5vs7)(co2 effect)', 'c.arabica DE GO (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora UP-DE VS c.arabica UP-DE (co2 effect):
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]),
                     set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']])],
                    ['c.canephora up-DE (5vs7)(co2 effect)', 'c.arabica up-DE (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora UP-DE GO VS c.arabica UP-DE GO (co2 effect):
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 5vs7, 1vs3 (GO > 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['5vs7']])),
                     set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['1vs3']]))],
                    ['c.canephora up-DE GO (5vs7)(co2 effect)', 'c.arabica up-DE GO (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora DOWN-DE VS c.arabica DOWN-DE (co2 effect):
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']])],
                    ['c.canephora down-DE (5vs7)(co2 effect)', 'c.arabica down-DE (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora DOWN-DE GO VS c.arabica DOWN-DE GO (co2 effect):
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 5vs7, 1vs3 (GO > 0)
        create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['5vs7']])),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['1vs3']]))],
                    ['c.canephora down-DE GO (5vs7)(co2 effect)', 'c.arabica down-DE GO (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora UP-DE VS c.canephora DOWN-DE:
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']])],
                    ['c.canephora UP-DE (5vs7)(co2 effect)', 'c.canephora DOWN-DE (5vs7)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora UP-DE GO VS c.canephora DOWN-DE GO:
        # results.csv 5vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 5vs7 (GO > 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['5vs7']])),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['5vs7']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['5vs7']]))],
                    ['c.canephora UP-DE GO (5vs7)(co2 effect)', 'c.canephora DOWN-DE GO (5vs7)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.arabica UP-DE VS c.arabica DOWN-DE:
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']]),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']])],
                    ['c.arabica UP-DE (1vs3)(co2 effect)', 'c.arabica DOWN-DE (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.arabica UP-DE GO VS c.arabica DOWN-DE GO:
        # results.csv 1vs3 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 1vs3 (GO > 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['1vs3']])),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs3']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['1vs3']]))],
                    ['c.arabica UP-DE GO (1vs3)(co2 effect)', 'c.arabica DOWN-DE GO (1vs3)(co2 effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 VS C700:
        # htseq_counts 1, 5 (count > 0)
        # htseq_counts 3, 7 (count > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 5]]),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3, 7]])],
                    ['C380 expressed genes (1and5)', 'C700 expressed genes (3and7)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 GO VS C700 GO:
        # htseq_counts 1, 5 (count > 0)
        # htseq_counts 3, 7 (count > 0)
        # GO ALL (GO > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1, 5]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()])),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3, 7]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()]))],
                    ['C380 expressed genes GO (1and5)', 'C700 expressed genes GO (3and7)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 DE VS C700 DE (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]),
                     set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']])],
                    ['C380 DE (1vs5)(genome effect)', 'C700 DE (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 DE GO VS C700 DE GO (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 5vs7, 1vs3 (GO > 0)
        create_venn([set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items() if key in ['1vs5']])),
                     set.union(*[_dict['up'].union(_dict['down']) for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items() if key in ['3vs7']]))],
                    ['C380 DE GO (1vs5)(genome effect)', 'C700 DE GO (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 UP-DE VS C700 UP-DE (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]),
                     set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']])],
                    ['C380 up-DE (1vs5)(genome effect)', 'C700 up-DE (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 UP-DE GO VS C700 UP-DE GO (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['1vs5']])),
                     set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['3vs7']]))],
                    ['C380 up-DE GO (1vs5)(genome effect)', 'C700 up-DE GO (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 DOWN-DE VS C700 DOWN-DE (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']])],
                    ['C380 down-DE (1vs5)(genome effect)', 'C700 down-DE (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 DOWN-DE GO VS C700 DOWN-DE GO (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['1vs5']])),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['3vs7']]))],
                    ['C380 down-DE GO (1vs5)(genome effect)', 'C700 down-DE GO (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 UP-DE VS C380 DOWN-DE (genome effect):
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']])],
                    ['C380 up-DE (1vs5)(genome effect)', 'C380 down-DE (1vs5)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 UP-DE GO VS C380 DOWN-DE GO:
        # results.csv 1vs5 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 1vs5 (GO > 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['1vs5']])),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['1vs5']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['1vs5']]))],
                    ['C380 up-DE GO (1vs5)(genome effect)', 'C380 down-DE GO (1vs5)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C700 UP-DE VS C700 DOWN-DE:
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']]),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']])],
                    ['C700 up-DE (3vs7)(genome effect)', 'C700 down-DE (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C700 UP-DE GO VS C700 DOWN-DE GO:
        # results.csv 3vs7 (FDR p-adjusted < 0.05, log2FoldChange > 0 & < 0)
        # GO 3vs7 (GO > 0)
        create_venn([set.union(*[_dict['up'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']]).intersection(set.union(*[_dict['up'] for key, _dict in seqtable_set_dict.items() if key in ['3vs7']])),
                     set.union(*[_dict['down'] for key, _dict in de_results_set_dict.items()
                                 if key in ['3vs7']]).intersection(set.union(*[_dict['down'] for key, _dict in seqtable_set_dict.items() if key in ['3vs7']]))],
                    ['C700 up-DE GO (3vs7)(genome effect)', 'C700 down-DE GO (3vs7)(genome effect)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.arabica C380 VS c.arabia C700:
        # htseq_counts 1, 3 (count > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1]]),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3]])],
                    ['c.arabica C380 expressed genes (1)', 'c.arabica C700 expressed genes (3)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.arabica C380 GO VS c.arabia C700 GO:
        # htseq_counts 1, 3 (count > 0)
        # GO ALL (GO > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()])),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()]))],
                    ['c.arabica C380 expressed genes GO (1)', 'c.arabica C700 expressed genes GO (3)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora C380 VS c.canephora C700:
        # htseq_counts 5, 7 (count > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [5]]),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [7]])],
                    ['c.canephora C380 expressed genes (5)', 'c.canephora C700 expressed genes (7)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # c.canephora C380 GO VS c.canephora C700 GO:
        # htseq_counts 5, 7 (count > 0)
        # GO ALL (GO > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [5]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()])),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [7]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()]))],
                    ['c.canephora C380 expressed genes GO (5)', 'c.canephora C700 expressed genes GO (7)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 c.arabica VS C380 c.canephora:
        # htseq_counts 1, 5 (count > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1]]),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [5]])],
                    ['C380 c.arabica expressed genes (1)', 'C380 c.canephora expressed genes (5)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C380 c.arabica GO VS C380 c.canephora GO:
        # htseq_counts 1, 5 (count > 0)
        # GO ALL (GO > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [1]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()])),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [5]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()]))],
                    ['C380 c.arabica expressed genes GO (1)', 'C380 c.canephora expressed genes GO (5)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C700 c.arabica VS C700 c.canephora:
        # htseq_counts 3, 7 (count > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3]]),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [7]])],
                    ['C700 c.arabica expressed genes (3)', 'C700 c.canephora expressed genes (7)'],
                    parser.output_folder,
                    ['png', 'eps'])

        # C700 c.arabica GO VS C700 c.canephora GO:
        # htseq_counts 3, 7 (count > 0)
        # GO ALL (GO > 0)
        create_venn([set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [3]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()])),
                     set.union(*[_set for key, _set in htseq_counts_set_dict.items() if key in [7]]).intersection(set.union(*[_dict['up'].union(_dict['down']) for key, _dict in seqtable_set_dict.items()]))],
                    ['C700 c.arabica expressed genes GO (3)', 'C700 c.canephora expressed genes GO (7)'],
                    parser.output_folder,
                    ['png', 'eps'])

        print('OK')

        #
        # values.csv
        #

        print('Generating values.csv ... ', end='')
        with open(os.path.join(parser.output_folder, 'values.csv'), 'w', newline='') as values_csv:
            values_csv_writer = csv.writer(values_csv)
            values_csv_writer.writerows(venn_values)
        print('OK')

    else:
        print('Please input valid and existing paths for \'htseq_counts_folder\',',
              '\'de_results_folder\', \'seqtable_folder\' and \'output_folder\'.')
