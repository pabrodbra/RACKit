from __future__ import print_function, division

taxa_level_dict = {
    "Specie": 6,
    "Genus": 5,
    "Family": 4,
    "Class" : 3,
    "Phylum": 2,
    "Order": 1,
    "Domain": 0
}

class TaxaLoad(object):
    def __init__(self, taxon_path_file, rank="Specie"):
        self.taxon_path = taxon_path_file
        self.taxa_rank = rank
        self.taxonomy_dictionary = {}
        self.create_tax_dictionary(rank)

    def create_tax_dictionary(self, rank="Specie"):
        self.taxonomy_dictionary = {}

        with open(self.taxon_path, 'r') as f:
            for line in f:
                items = line.replace('\n', '').split[',']
                seq_id = items[0]
                try:
                    assigned_taxa = items[1].replace('"', '').split(';')[taxa_level_dict[rank]]
                    self.taxonomy_dictionary[seq_id] = assigned_taxa
                except:
                    pass

    def save_taxonomy_dictionary(self, out_name="ID-Category.csv"):
        with open(out_name, 'w') as f:
            # for key in sorted(self.taxonomy_dictionary):
            for key in self.taxonomy_dictionary.keys():
                current_tax = self.taxonomy_dictionary.get(key)
                f.write(str(key) + "," + str(current_tax) + "\n")
    
    def get_taxonomy_dictionary(self):
        return self.taxonomy_dictionary