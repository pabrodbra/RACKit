from __future__ import print_function
import sys


"""
GREFCO - Gather REads Forming Contigs
Developed by: Pablo Rodriguez Brazzarola
"""


class GREFCO(object):
    def __init__(self, file_path):
        self.grefco_file = file_path
        self.parsed_file = file_path
        self.contig_dictionary = {}
        self.reads_dictionary = {}

    def make_full_contig_dictionary(self):
        if self.parsed_file is "":
            print("No Filter File specified")
            exit(0)

        with open(self.parsed_file, 'r') as f:
            for line in f:
                line = line.replace("\n", "").replace(">", "").replace(" ", ";")
                split = line.split(";")
                contig_id = split[1]
                split.pop(1)
                values = split

                if contig_id not in self.contig_dictionary:
                    reads = list()
                    reads.append(values)
                    self.contig_dictionary[contig_id] = reads
                else:
                    reads = self.contig_dictionary.get(contig_id)
                    reads.append(values)

    def make_full_read_dictionary(self):
        if self.parsed_file is "":
            print("No Filter File specified")
            exit(0)

        with open(self.parsed_file, 'r') as f:
            for line in f:
                line = line.replace("\n", "").replace(">", "").replace(" ", ";")
                split = line.split(";")
                read_id = split[0]
                split.pop(0)
                values = split

                if read_id not in self.reads_dictionary:
                    contigs = list()
                    contigs.append(values)
                    self.reads_dictionary[read_id] = contigs
                else:
                    contigs = self.reads_dictionary.get(read_id)
                    contigs.append(values)

    @staticmethod
    def output_dictionary(output, dictionary, id_flag=False):
        if id_flag is True:
            index_f = open(output + ".index", 'w')
        with open(output, 'w') as f:
            key_counter = 0
            for key in dictionary.keys():
                seq_id = key
                if id_flag is True:
                    index_f.write(str(seq_id) + "," + str(key_counter) + '\n')
                    seq_id = str(seq_id) + "_" + str(key_counter)
                    # seq_id = str(key_counter) + ";" + seq_id
                values = dictionary.get(key)
                line = str(seq_id) + ";" + str(len(values)) + ";"
                for val in values:
                    line += str(val) + ","
                line = line[:-1] + "\n"
                f.write(line)
                key_counter += 1

        if id_flag is True:
            index_f.close()
