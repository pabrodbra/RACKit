from __future__ import print_function, division
import glob


class TaxonomyLoad(object):
    def __init__(self, folder):
        self.taxonomy_folder = folder
        self.taxonomy_dictionary = {}
        self.index_dictionary = {}

    def create_tax_dictionary(self):
        self.taxonomy_dictionary = {}

        path = self.taxonomy_folder + "*.fasta"
        files = glob.glob(pathname=path)

        for tax_file in files:
            # WINDOWS: tax_name = tax_file.split('\\')[-1].replace("SUMreads-", "").replace(".fasta", "")
            # UNIX:
            tax_name = tax_file.split('/')[-1].replace("SUMreads-", "").replace(".fasta", "")
            with open(tax_file, 'r') as f:
                for line in f:
                    seq_id = line.replace(">", "").replace("\n", "")
                    self.taxonomy_dictionary[seq_id] = tax_name

    def create_tax_dictionary_indexed(self):
        self.taxonomy_dictionary = {}
        f_index = open("Category-Index.csv", 'w')

        path = self.taxonomy_folder + "*.fasta"
        files = glob.glob(pathname=path)

        tax_index = 0
        for tax_file in files:
            tax_name = tax_file.split('\\')[-1].replace("SUMreads-", "").replace(".fasta", "")
            f_index.write(tax_name + "," + str(tax_index) + "\n")
            with open(tax_file, 'r') as f:
                for line in f:
                    seq_id = line.replace(">", "").replace("\n", "")
                    self.taxonomy_dictionary[seq_id] = tax_index
            tax_index += 1

        f_index.close()

    def save_taxonomy_dictionary(self, out_name="ID-Category.csv"):
        with open(out_name, 'w') as f:
            # for key in sorted(self.taxonomy_dictionary):
            for key in self.taxonomy_dictionary.keys():
                current_tax = self.taxonomy_dictionary.get(key)
                f.write(str(key) + "," + str(current_tax) + "\n")

    def get_taxonomy(self, id):
        if id in self.get_taxonomy_dictionary().keys():
            return self.get_taxonomy_dictionary().get(id)
        else:
            return None

    @staticmethod
    def sum_indexes(index1, index2, out_file="SUM-indexes.csv"):
        f1 = open(index1, 'r')
        f2 = open(index2, 'r')
        with open(out_file, 'w') as f:
            for line in f1:
                f.write(line)
            for line in f2:
                f.write(line)
        f1.close()
        f2.close()

    def load_indexes(self, index_name="SUM-indexes.csv"):
        num_lines = sum(1 for line in open(index_name, 'r'))
        plus_line_counter = 1
        with open(index_name, 'r') as f:
            for line in f:
                split = line.split(',')
                key_id = split[0]
                key_index = split[1]
                indexed_key = str(key_id) + "_" + str(key_index)
                tax_is = self.get_taxonomy(key_id)
                self.get_taxonomy_dictionary().pop(key_id, None)
                if tax_is is not None:
                    self.get_taxonomy_dictionary()[indexed_key] = tax_is
                else:
                    indexed_key = str(key_id) + "_" + str(num_lines + plus_line_counter)
                    self.get_taxonomy_dictionary()[indexed_key] = tax_is
                    plus_line_counter += 1

    def load_taxonomy_dictionary(self, in_name="ID-Category.csv"):
        self.taxonomy_dictionary = {}

        with open(in_name, 'r') as f:
            for line in f:
                split = line.split(',')
                self.taxonomy_dictionary[split[0]] = split[1]

    def get_taxonomy_dictionary(self):
        return self.taxonomy_dictionary

if __name__ == "__main__":
    read_taxonomy = TaxonomyLoad(
        # "3 - Family/Reads by Family/")
        "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/Reads by Family/")
    contig_taxonomy = TaxonomyLoad(
        # "3 - Family/Contigs by Family/")
        "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/Contigs by Family/")

    #read_taxonomy.create_tax_dictionary()
    #contig_taxonomy.create_tax_dictionary()

    read_taxonomy.create_tax_dictionary_indexed()
    contig_taxonomy.create_tax_dictionary_indexed()

    #read_taxonomy.load_taxonomy_dictionary()
    #contig_taxonomy.load_taxonomy_dictionary()

    print("Taxons found for reads:")
    print(len(read_taxonomy.get_taxonomy_dictionary().keys()))

    print("Taxons found for contigs:")
    print(len(contig_taxonomy.get_taxonomy_dictionary().keys()))

    read_taxonomy.save_taxonomy_dictionary(out_name="ReadID-CategoryIndex.csv")
    contig_taxonomy.save_taxonomy_dictionary(out_name="ContigID-CategoryIndex.csv")

    """
    forward_reads_grefco = GREFCO("data/forward_parsed-reads.grefco")
    reverse_reads_grefco = GREFCO("data/reverse_parsed-reads.grefco")

    forward_reads_grefco.load_grefco(type="reads")
    reverse_reads_grefco.load_grefco(type="reads")
    ---------------
    index_1 = "data/forward_parsed-reads-indexedplus.grefco.index"
    index_2 = "data/reverse_parsed-reads-indexedplus.grefco.index"
    index_summed = "data/SUM-indexes.csv"
    # TaxonomyLoad.sum_indexes(index1=index_1, index2=index_2, out_file="data/SUM-indexes.csv")

    # read_taxonomy.load_indexes(index_summed)
    """