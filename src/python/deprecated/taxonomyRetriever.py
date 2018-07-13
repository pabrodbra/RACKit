from __future__ import print_function, division
import sys
import glob

from grefco import GREFCO


class TaxonomyRetriever(object):
    def __init__(self, grefco_dictionary):
        self.sequence_dictionary = grefco_dictionary

    def get_all_taxonomies(self, folders, name="Taxonomy.txt"):
        path = folders + "*.fasta"
        files = glob.glob(pathname=path)
        all_taxonomies = []
        total_taxonomies = 0

        for f in files:
            temp_res = self.get_taxonomy_count(f)
            new_res = (temp_res[0].split('\\')[-1], temp_res[1])
            all_taxonomies.append(new_res)
            total_taxonomies += new_res[1]

        TaxonomyRetriever.create_taxonomy_files(all_taxonomies, name)

        return all_taxonomies, total_taxonomies

    def get_taxonomy_count(self, taxonomy_file):
        raise NotImplementedError("Please implement this method")

    @staticmethod
    def create_taxonomy_files(list_taxonomies, name):
        with open(name, 'w') as f:
            for tax_info in list_taxonomies:
                f.write('"' + tax_info[0] + '",' + str(tax_info[1]) + "\n")

        return name

    @staticmethod
    def sum_taxonomy_files(list_files, name="SumTax.csv"):
        dict_taxons = {}
        tot_len = 0

        for tax_file in list_files:
            with open(tax_file, 'r') as f:
                for line in f:
                    split = line.replace("\n", "").split(",")

                    if dict_taxons.has_key(split[0]):
                        dict_taxons[split[0]] += int(split[1])
                    else:
                        dict_taxons[split[0]] = int(split[1])

        with open(name, 'w') as f:
            for key in sorted(dict_taxons):
                f.write(key + ',' + str(dict_taxons[key]) + "\n")
                tot_len += int(dict_taxons[key])

        return tot_len

    @staticmethod
    def normalize_taxonomy_files(total_sequences, name="SumTax.csv"):
        with open(name + "-normalized", 'w') as out_f:
            with open(name, 'r') as f:
                for line in f:
                    split = line.replace("\n", "").split(",")
                    norm_val = int(split[1])/total_sequences * 100
                    out_f.write(split[0] + "," + str(norm_val) + "\n")


class SimpleContigTaxonomyRetriever(TaxonomyRetriever):
    def get_taxonomy_count(self, taxonomy_file):
        """
        Retrieve Taxonomy Name  + Amount of contigs associated with a contigs taxonomy.
        GREFCO must have had been declared as 'simple'
        :param taxonomy_file: File that lists all contigs associated with a contigs taxonomy file
        :return: Tuple as (Taxonomy Name, Taxonomy Count)
        """
        taxonomy_name = taxonomy_file.replace("SUMreads-", "").replace(".fasta", "")
        taxonomy_count = 0
        with open(taxonomy_file, 'r') as f:
            for line in f:
                line = line.replace("\n", "").replace(">", "")
                if self.sequence_dictionary.has_key(line):
                    temp_val = self.sequence_dictionary.get(line)
                    taxonomy_count = taxonomy_count + len(temp_val)

        taxonomy_name = taxonomy_name.split('/')[-1]
        taxonomy_info = (taxonomy_name, taxonomy_count)
        return taxonomy_info


# FIX
class ComplexContigTaxonomyRetriever(TaxonomyRetriever):
    def get_taxonomy_count(self, taxonomy_file):
        """
        Retrieve Taxonomy Name + Amount of contigs associated with a contigs taxonomy.
        GREFCO must have had been declared as 'complex'
        :param taxonomy_file: File that lists all contigs associated with that contigs taxonomy file
        :return: Tuple as (Taxonomy Name, Taxonomy Count)
        """
        taxonomy_name = taxonomy_file.replace("SUMreads-", "").replace(".fasta", "")
        taxonomy_count = 0
        print(taxonomy_file)
        with open(taxonomy_file, 'r') as f:
            for line in f:
                seq_id = line.replace("\n", "").replace(">", "")
                if self.sequence_dictionary.has_key(line):
                    temp_val = self.sequence_dictionary.get(seq_id)
                    best_score = 1.0
                    for item in temp_val:

                        if best_score > item:
                            best_score = item[1]
                    taxonomy_count += 1

        taxonomy_name = taxonomy_name.split('/')[-1]
        taxonomy_info = (taxonomy_name, taxonomy_count)
        return taxonomy_info


# FIX
class SimpleReadTaxonomyRetriever(TaxonomyRetriever):
    def get_taxonomy_count(self, taxonomy_file):
        """
        Retrieve Taxonomy Name  + Amount of reads associated with a contigs taxonomy.
        GREFCO must have had been declared as 'simple'
        :param taxonomy_file: File that lists all contigs associated with that taxonomy file
        :return: Tuple as (Taxonomy Name, Taxonomy Count)
        """
        taxonomy_name = taxonomy_file.replace("SUMreads-", "").replace(".fasta", "")
        taxonomy_count = 0
        with open(taxonomy_file, 'r') as f:
            for line in f:
                line = line.replace("\n", "").replace(">", "")
                if self.sequence_dictionary.has_key(line):
                    temp_val = self.sequence_dictionary.get(line)
                    taxonomy_count = taxonomy_count + len(temp_val)

        taxonomy_name = taxonomy_name.split('/')[-1]
        taxonomy_info = (taxonomy_name, taxonomy_count)
        return taxonomy_info


# FIX
class ComplexReadTaxonomyRetriever(TaxonomyRetriever):
    def get_taxonomy_count(self, taxonomy_file):
        """
        Retrieve Taxonomy Name + Amount of contigs associated with a contigs taxonomy.
        GREFCO must have had been declared as 'complex'
        :param taxonomy_file: File that lists all contigs associated with that contigs taxonomy file
        :return: Tuple as (Taxonomy Name, Taxonomy Count)
        """
        taxonomy_name = taxonomy_file.replace("SUMreads-", "").replace(".fasta", "")
        taxonomy_count = 0
        with open(taxonomy_file, 'r') as f:
            for line in f:
                seq_id = line.replace("\n", "").replace(">", "")
                values = self.sequence_dictionary.keys()
                print(values)
                for value in values:
                    if value[0] == seq_id:
                        temp_val = self.sequence_dictionary.get(k)
                        print(temp_val)
                        best_score = 1.0
                        for item in temp_val:
                            if best_score > item[1]:
                                best_score = item[1]
                        taxonomy_count += 1

        taxonomy_name = taxonomy_name.split('/')[-1]
        taxonomy_info = (taxonomy_name, taxonomy_count)
        return taxonomy_info


if __name__ == "__main__":
    folder = "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/Contigs by Family/"
    '''    '''
    # ########### TEST SUM AND NORMALIZE #############
    tax_names = ["data/forward.filtered.grefco-Taxonomy.csv", "data/reverse.filtered.grefco-Taxonomy.csv"]

    sum_tax_name = "data/SumTax.csv"
    total_seqs = TaxonomyRetriever.sum_taxonomy_files(tax_names, sum_tax_name)
    print(total_seqs)
    TaxonomyRetriever.normalize_taxonomy_files(total_sequences=total_seqs, name=sum_tax_name)

    ''' 
    if len(sys.argv) < 2:
        print("USAGE: python taxonomyRetriever.py <grefco file>", file=sys.stderr)
        exit(-1)
    '''

    '''
    # ########### TEST LOAD AND CREATE #############
    
    grefco1 = GREFCO(sys.argv[1])
    # grefco2 = GREFCO(sys.argv[2])

    grefco1.read_grefco()
    # grefco2.read_grefco()

    taxonomy_retriever1 = SimpleContigTaxonomyRetriever(grefco1.contig_dictionary)
    # taxonomy_retriever1 = ComplexReadTaxonomyRetriever(grefco1.reads_dictionary)
    # taxonomy_retriever2 = TaxonomyRetriever(grefco2)

    test1, test_total1 = taxonomy_retriever1.get_all_taxonomies(folder, sys.argv[1] + "-Taxonomy.csv")
    # test2, test_total2 = taxonomy_retriever2.get_all_taxonomies(folder, sys.argv[1] + "-Taxonomy.csv")

    print("Test Total Length 1 = " + str(test_total1))
    # print("Test Total Length 2 = " + str(test_total2)) '''


# python taxonomyRetriever.py data/forward.filtered.grefco
'''

    
'''