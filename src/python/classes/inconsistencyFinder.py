from __future__ import print_function, division
import concurrent.futures

from classes.grefco import GREFCO
from classes.taxonomyLoad import TaxonomyLoad


class InconsistencyFinder(object):
    def __init__(self, reads_dict, reads_tax, contigs_tax):
        self.reads_dictionary = reads_dict
        self.reads_taxonomy = reads_tax
        self.contigs_taxonomy = contigs_tax

    def get_read_tax(self, read_id):
        ret = None
        if read_id in self.reads_taxonomy.keys():
            ret = self.reads_taxonomy.get(read_id)
        return ret

    def get_contig_tax(self, contig_id):
        ret = None
        if contig_id in self.contigs_taxonomy.keys():
            ret = self.contigs_taxonomy.get(contig_id)
        return ret

    def write_inconsistencies_sequential(self, inconsistency_fname="REVCO.sinconsistency"):
        inconsistency_counter = 0
        temp_sc = ''

        with open(inconsistency_fname, 'w') as f:
            f.write("ReadID,ReadTaxonomy,ContigID,ContigTaxonomy,SemiConsistant\n")
            for read_id in self.reads_dictionary.keys():
                new_inconsistency, temp_sc = self.find_inconsistency_sequential(read_id)

                if len(new_inconsistency) == 0:
                    pass
                else:
                    for ni in new_inconsistency:
                        f.write(str(ni[0]) + "," + str(ni[1])
                                + "," + str(ni[2]) + "," + str(ni[3]) + "," + temp_sc + "\n")
                        inconsistency_counter += 1

        return inconsistency_counter

    def find_inconsistency_sequential(self, read_id):
        inconsistencies_found = []
        semi_consistant = '-'

        read_taxon = self.get_read_tax(read_id)
        contig_values = self.reads_dictionary.get(read_id)

        for contig_match in contig_values:
            contig_id = contig_match[0]
            contig_tax = self.get_contig_tax(contig_id)
            if read_taxon != contig_tax:
                new_info = (str(read_id), str(read_taxon), str(contig_id), str(contig_tax))
                inconsistencies_found.append(new_info)
            else:
                semi_consistant = '+' # At least one match correct

        return inconsistencies_found, semi_consistant