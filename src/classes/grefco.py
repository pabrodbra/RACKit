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

    # FIX
    def load_full_dictioanry(self, fname, type='contig'):
        new_dict = dict()

        with open(fname, 'r') as f:
            lines = f.readlines()
            for line in lines:
                items = line.rsplit(';')
                id = items[0]
                values = items[2].replace('[', '').replace('], ', ";").replace(']', '').rsplit(';')

                new_values = list()
                for value in values:
                    new_items = value.rsplit(',')
                    new_values.append(new_items)

                new_dict[id] = new_values

        if type == 'contig':
            self.contig_dictionary = new_dict
        elif type == 'read':
            self.reads_dictionary = new_dict

    # OLD

    def make_simple_dictionary(self):
        if self.parsed_file is "":
            print("No Parsed File specified")
            exit(0)

        with open(self.parsed_file, 'r') as f:
            for line in f:
                # line = line.replace("\n", "").replace(">", "")
                split = line.split(" ")
                if not self.contig_dictionary.has_key(split[0]):
                    reads = list()
                    reads.append(split[1])
                    self.contig_dictionary[split[0]] = reads
                elif self.contig_dictionary.has_key(split[0]):
                    reads = self.contig_dictionary.get(split[0])
                    reads.append(split[1])

    def filter_for_complex_dictionary(self):
        if self.parsed_file is "":
            print("No Filter File specified")
            exit(0)

        filtered_lines = []
        with open(self.parsed_file, 'r') as f:
            lines = ""
            open_set = False
            for line in f:
                if line[0] == ">":
                    if open_set is True:
                        filtered_lines.append(lines)
                        lines = ""
                        open_set = False
                    if open_set is False:
                        split = line.split(" ")
                        new_line = split[2] + " " + split[0]
                        lines += new_line
                        open_set = True
                elif line[0] != ">" and open_set is True:
                    split = line.split("\t")
                    new_line = " " + split[2] + " " + split[3] + " " + split[4]
                    lines += new_line

        return filtered_lines

    def make_complex_dictionary(self):
        filtered_lines = self.filter_for_complex_dictionary()
        for line in filtered_lines:
            line = line.replace("\n", "").replace(">", "")
            split = line.split(" ")
            key = split[0]

            counter_next_hit = 2
            best_triple = (1, 1, 1)
            max_info = (len(split)-2)/3

            for i in range(0,max_info,1):
                current_triple = (split[counter_next_hit], split[counter_next_hit + 1], split[counter_next_hit + 2])
                if best_triple[0] < current_triple:
                    best_triple = current_triple
                counter_next_hit += 3

            values = (split[1], best_triple[0], best_triple[1], best_triple[2])
            if not self.contig_dictionary.has_key(key):
                reads = list()
                reads.append(values)
                self.contig_dictionary[split[0]] = reads
            elif self.contig_dictionary.has_key(key):
                reads = self.contig_dictionary.get(key)
                reads.append(values)

    def load_grefco(self, type="contig"):
        if self.grefco_file is "":
            print("No GREFCO File specified")
            exit(0)

        with open(self.grefco_file, 'r') as f:
            for line in f:
                line = line.replace("\n", "")
                split = line.split(";")
                if type == "contigs":
                    self.contig_dictionary[split[0]] = split[2].split(",")
                elif type == "reads":
                    self.reads_dictionary[split[0]] = split[2].split(",")

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

    def invert_simple_dictionary(self):
        self.reads_dictionary = dict()
        for k, values in self.contig_dictionary.iteritems():
            for v in values:
                if not self.reads_dictionary.has_key(v):
                    new_val = list()
                    new_val.append(k)
                    self.reads_dictionary[v] = new_val
                elif self.reads_dictionary.has_key(v):
                    reads = self.reads_dictionary.get(v)
                    reads.append(k)

    def invert_complex_dictionary(self):
        self.reads_dictionary = dict()
        for k, values in self.contig_dictionary.iteritems():
            for v in values:
                if not self.reads_dictionary.has_key(v[0]):
                    new_val = list()
                    new_val_info = (k, v[1], v[2], v[3])
                    new_val.append(new_val_info)
                    self.reads_dictionary[v[0]] = new_val
                elif self.reads_dictionary.has_key(v[0]):
                    reads = self.reads_dictionary.get(v[0])
                    new_val_info = (k, v[1], v[2], v[3])
                    reads.append(new_val_info)

# python grefco.py data/forward.filtered

if __name__ == "__main__":
    ''' 
    if len(sys.argv) != 2:
        print("USAGE: python grefco.py <parsed blast file>", file=sys.stderr)
        exit(-1)
    '''
    # name = sys.argv[1]

    # name = "data/forward.filtered"
    # name = "data/reverse.filtered"
    # grefco.make_simple_dictionary()

    name = "data/forward_parsed"
    name2 = "data/reverse_parsed"

    grefco = GREFCO(name)
    grefco2 = GREFCO(name2)
    grefco.make_complex_dictionary()
    grefco2.make_complex_dictionary()

    output_name = name + ".grefco"
    output_name2 = name2 + ".grefco"

    GREFCO.output_dictionary(output_name, grefco.contig_dictionary)
    GREFCO.output_dictionary(output_name2, grefco2.contig_dictionary)

    # refco.invert_simple_dictionary()
    grefco.invert_complex_dictionary()
    grefco2.invert_complex_dictionary()
    GREFCO.output_dictionary(name + "-reads-indexedplus.grefco", grefco.reads_dictionary, id_flag=True)
    GREFCO.output_dictionary(name2 + "-reads-indexedplus.grefco", grefco2.reads_dictionary, id_flag=True)
    # GREFCO.output_dictionary(sys.argv[1] + "-reads-indexed.grefco", grefco.reads_dictionary)

'''
PLOT (ARGUMENTOS)
INCONSISTENCIA
PLOT (FROM COMPLEX READ DICTIONARY -> GET MATCH WHERE READ HITS WITH HIGHEST SCORE -> GET TAX OF MATCH -> SUM TO TAXFILE
'''

'''
Simple Dictionary : Requires filtered + parsed
Complex Dictionary : Requires parsed
------
Indexed : For C++ InconsistencyFinder and Add TaxonomyGroup in C++ per index
    If no TaxonomyGroup, set to "NA"
When SAVE INDEXED PLUS DICTIONARY , SAVE HIGHEST INDEX for allocating memory later
'''