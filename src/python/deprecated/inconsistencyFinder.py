from __future__ import print_function, division
import concurrent.futures

from grefco import GREFCO
from taxonomyLoad import TaxonomyLoad


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

                # Sequential
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
                semi_consistant = '+'

        return inconsistencies_found, semi_consistant

    def write_inconsistencies_parallel(self, inconsistency_fname="REVCO.pinconsistency", num_workers=8, m_partitions=5):
        inconsistency_counter = 0

        with open(inconsistency_fname, 'w') as f:
            f.write("ReadID,ReadTaxonomy,ContigID,ContigTaxonomy\n")

            read_blocks = self.split_reads(num_workers, m_partitions)

            with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                future_to_readblock = {executor.submit(self.find_inconsistency_parallel, read_block):
                                           read_block for read_block in read_blocks}
                for future in concurrent.futures.as_completed(future_to_readblock):
                    rb_id = future_to_readblock[future]
                    new_inconsistencies = future.result()

                    for ni in new_inconsistencies:
                        if len(ni) == 0:
                            pass
                        else:
                            f.write(str(ni[0]) + "," + str(ni[1])
                                    + "," + str(ni[2]) + "," + str(ni[3]) + "\n")
                            inconsistency_counter += 1

        return inconsistency_counter

    def split_reads(self, num_threads=8, n_partitions=4):
        read_keys = list(self.reads_dictionary.keys())
        blocks = list()
        n_reads = len(read_keys)

        m_partitions = int(n_reads/n_partitions)
        even_chunks = [read_keys[i:i+m_partitions] for i in range(0, n_reads, m_partitions)]
        current_split = 1

        for chunk in even_chunks:
            num_partitions = int(n_reads/(num_threads*current_split))
            current_blocks = [chunk[i:i+num_partitions] for i in range(0, len(chunk), num_partitions)]
            blocks.extend(current_blocks)
            current_split += 1

        """ split_begin, split_reduce_factor, 
        reads_in_blocks = 0
        previous_block_limit = 0
        split_counter = 0
        max_split = split_reduce_factor * 8
        current_split = split_begin * 2

        print("Begin Split")
        current_block_size = int(n_reads / current_split)
        while reads_in_blocks < n_reads:
            new_block_limit = previous_block_limit + current_block_size
            new_block = read_keys[previous_block_limit:new_block_limit]

            blocks.append(new_block)

            previous_block_limit = reads_in_blocks = new_block_limit
            split_counter += 1

            '''
            if split_counter < split_reduce_factor and current_split < max_split:
                current_block_size = int(n_reads / current_split)
                new_block_limit = previous_block_limit + current_block_size
                new_block = read_keys[previous_block_limit:new_block_limit]

                blocks.append(new_block)

                previous_block_limit = reads_in_blocks = new_block_limit
                split_counter += 1
            elif split_counter >= split_reduce_factor:
                if current_split < max_split:
                    current_split = current_split * 2
                split_counter = 0
                #print("--- Reducing Split: " + str(current_split))
            '''

        print("Finished Split")
        """
        return blocks

    def find_inconsistency_parallel(self, read_block):
        inconsistencies_found = []

        for read_id in read_block:
            read_taxon = self.get_read_tax(read_id)
            contig_values = self.reads_dictionary.get(read_id)

            for contig_match in contig_values:
                contig_id = contig_match[0]
                contig_tax = self.get_contig_tax(contig_id)
                if read_taxon != contig_tax:
                    new_info = (str(read_id), str(read_taxon), str(contig_id), str(contig_tax))
                    inconsistencies_found.append(new_info)

        return inconsistencies_found

    # OLD
    def write_inconsistencies(self, inconsistency_fname="REVCO.inconsistency"):
        inconsistency_counter = 0

        with open(inconsistency_fname, 'w') as f:
            f.write("ReadID,ReadTaxonomy,ContigID,ContigTaxonomy\n")
            for read_id in self.reads_dictionary.keys():
                read_tax = self.get_read_tax(read_id)
                values = self.reads_dictionary.get(read_id)
                current_id_index = 0

                # if read_tax is not None and (len(values) / 4) > 1:
                if (len(values) / 4) > 1:
                    while current_id_index < len(values):
                        value = values[current_id_index].replace("(", "").replace("'", "")
                        contig_id = value
                        contig_tax = self.get_contig_tax(contig_id)

                        # Output the inconsistency
                        if contig_tax is not None and read_tax != contig_tax:
                            f.write(str(read_id) + "," + str(read_tax)
                                    + "," + str(contig_id) + "," + str(contig_tax) + "\n")
                            inconsistency_counter += 1
                        current_id_index += 4
                # elif read_tax is not None and len(values) / 4 == 1:
                elif len(values) / 4 == 1:
                    value = values[current_id_index].replace("(", "").replace("'", "")
                    contig_id = value
                    contig_tax = self.get_contig_tax(contig_id)
                    if contig_tax is not None and read_tax != contig_tax:
                        # Output the inconsistency
                        f.write(str(read_id) + "," + str(read_tax)
                                + "," + str(contig_id) + "," + str(contig_tax) + "\n")
                        inconsistency_counter += 1
                    current_id_index += 4

        return inconsistency_counter


# -------------------------

if __name__ == "__main__":
    forward_reads_grefco = GREFCO("data/forward_parsed-reads.grefco")
    reverse_reads_grefco = GREFCO("data/reverse_parsed-reads.grefco")

    forward_reads_grefco.load_grefco(type="reads")
    reverse_reads_grefco.load_grefco(type="reads")
    read_taxonomy = TaxonomyLoad(
        # "3 - Family/Reads by Family/")
        "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/Reads by Family/")
    contig_taxonomy = TaxonomyLoad(
        # "3 - Family/Contigs by Family/")
        "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/MEGAN/3 - Family/Contigs by Family/")
    """
    read_taxonomy.create_tax_dictionary()
    contig_taxonomy.create_tax_dictionary()

    #read_taxonomy.load_taxonomy_dictionary()
    #contig_taxonomy.load_taxonomy_dictionary()

    print("Taxons found for reads:")
    print(len(read_taxonomy.get_taxonomy_dictionary().keys()))

    print("Taxons found for contigs:")
    print(len(contig_taxonomy.get_taxonomy_dictionary().keys()))

    read_taxonomy.save_taxonomy_dictionary(out_name="ReadID-Category.csv")
    contig_taxonomy.save_taxonomy_dictionary(out_name="ContigID-Category.csv")
    """
    forward_inconsistency_finder = InconsistencyFinder(forward_reads_grefco.reads_dictionary,
                                                       read_taxonomy.get_taxonomy_dictionary(),
                                                       contig_taxonomy.get_taxonomy_dictionary())

    total_inconsistencies1 = forward_inconsistency_finder.write_inconsistencies(
        inconsistency_fname="forward-REVCO.inconsistency")

    reverse_inconsistency_finder = InconsistencyFinder(reverse_reads_grefco.reads_dictionary,
                                                       read_taxonomy.get_taxonomy_dictionary(),
                                                       contig_taxonomy.get_taxonomy_dictionary())

    total_inconsistencies2 = reverse_inconsistency_finder.write_inconsistencies(
        inconsistency_fname="reverse-REVCO.inconsistency")

    with open("total_inconsistency.txt", 'w') as f:
        f.write("Total number of Inconsistencies (Forward) Found: " + str(total_inconsistencies1) + "\n")
        f.write("Total number of Inconsistencies (Reverse) Found: " + str(total_inconsistencies2))
    print("Total number of Inconsistencies (Forward) Found: " + str(total_inconsistencies1))

    print("Total number of Inconsistencies (Reverse) Found: " + str(total_inconsistencies2))
"""
"""

'''
Requires: Folder from MEGAN Summarized IDs of Reads/Contigs for all Taxonomies

'''
