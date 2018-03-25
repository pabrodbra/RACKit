from __future__ import print_function, division

taxa_level_dict = {
    "Specie": -1,
    "Genus": -2,
    "Family": -3,
    "Class" : -4,
    "Phylum": -5,
    "Order": -6,
    "Domain": -7
}

class TaxaLoad(object):
    def __init__(self, taxon_path_file, rank="Specie"):
        self.taxon_path = taxon_path_file
        self.taxa_rank = rank
        self.taxonomy_dictionary = {}
        self.create_tax_dictionary(rank)
        self.inverted_dictionary={}

    def create_tax_dictionary(self, rank="Specie"):
        self.taxonomy_dictionary = {}

        with open(self.taxon_path, 'r') as f:
            for line in f:
                items = line.replace('\n', '').split(',')
                seq_id = items[0]
                try:
                    assigned_taxa = items[1].replace('"', '').split(';')[:-1][taxa_level_dict[rank]]
                    self.taxonomy_dictionary[seq_id] = assigned_taxa
                except:
                    pass

    def save_taxonomy_dictionary(self, out_name="ID-Category.csv"):
        with open(out_name, 'w') as f:
            # for key in sorted(self.taxonomy_dictionary):
            for key in self.taxonomy_dictionary.keys():
                current_tax = self.taxonomy_dictionary.get(key)
                f.write(str(key) + "," + str(current_tax) + "\n")
    
    def create_inverted_dictionary(self):
        self.inverted_dictionary = {v: k for k, v in self.taxonomy_dictionary.items()}

    def save_inverted_dictionary_count(self, output="Category-Count.csv"):
        with open(output, 'w') as f:
            for key in self.inverted_dictionary.keys():
                print(key); print(self.inverted_dictionary[key]);f.write( str(key) + ',' + str(len(self.inverted_dictionary[key])) + '\n' )

    def get_taxonomy_dictionary(self):
        return self.taxonomy_dictionary
