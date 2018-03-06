from classes.grefco import GREFCO
from classes.taxonomyLoad import TaxonomyLoad
from classes.inconsistencyFinder import InconsistencyFinder

import sys
import time

#import taxonomyRetriever

def grefco_main():
    filtered_file_name = sys.argv[1]
    #filtered_file_name = "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/test_parsed"
    #filtered_file_name = "/home/pablorod/data/metagenomics/REVCO/filter1_CONTIGBLAST_format_reads.pblast"
    grefco = GREFCO(filtered_file_name)

    grefco.make_full_read_dictionary()
    grefco.make_full_contig_dictionary()

    #grefco.load_full_dictioanry(filtered_file_name + ".grefcocontig", 'contig')
    #grefco.load_full_dictioanry(filtered_file_name + ".grefcoread", 'read')

    grefco_output = filtered_file_name + ".grefco"

    GREFCO.output_dictionary(grefco_output + "contig", grefco.contig_dictionary)

    GREFCO.output_dictionary(grefco_output + "read", grefco.reads_dictionary)

    return grefco


def taxonomy_load_main(tax_level="Species"):
    read_taxonomy = TaxonomyLoad(
        sys.argv[2] + tax_level + "/Reads/")
        #"/home/pablorod/data/metagenomics/REVCO/MEGAN/" + tax_level + "/Reads/")
        # "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/MEGAN/Species/ReadsT/")
    contig_taxonomy = TaxonomyLoad(
        sys.argv[2] + tax_level + "/Contigs/")
        #"/home/pablorod/data/metagenomics/REVCO/MEGAN/" + tax_level + "/Contigs/")
        # "C:/Users/Blinsky.Blinsk/Documents/- UMA/BitLab/Metagenome/Synthetic Dataset/MEGAN/Species/ContigNS/")

    read_taxonomy.create_tax_dictionary()
    contig_taxonomy.create_tax_dictionary()

    #read_taxonomy.create_tax_dictionary_indexed()
    #contig_taxonomy.create_tax_dictionary_indexed()

    #read_taxonomy.load_taxonomy_dictionary("ReadID-CategoryIndex.csv")
    #contig_taxonomy.load_taxonomy_dictionary("ContigID-CategoryIndex.csv")
    print("----- Comparing @ <" + tax_level + "> -----")
    print("Taxons found for reads:")
    print(len(read_taxonomy.get_taxonomy_dictionary().keys()))

    print("Taxons found for contigs:")
    print(len(contig_taxonomy.get_taxonomy_dictionary().keys()))

    read_taxonomy.save_taxonomy_dictionary(out_name="ReadID-"+tax_level+"Index.csv")
    contig_taxonomy.save_taxonomy_dictionary(out_name="ContigID-"+tax_level+"Index.csv")

    return read_taxonomy, contig_taxonomy


def inconsistency_finder_main(grefco, read_tax, contig_tax, out_name="inconsistency.csv"):
    inconsistency_finder = InconsistencyFinder(grefco.reads_dictionary,
                                                       read_tax.get_taxonomy_dictionary(),
                                                       contig_tax.get_taxonomy_dictionary())

    seq_i = time.time()
    total_inconsistencies = inconsistency_finder.write_inconsistencies_sequential(
        inconsistency_fname=out_name)
    seq_f = time.time()

    print("Total number of Inconsistencies Found: " + str(total_inconsistencies))
    print("Sequential Time: " + str(seq_f-seq_i))

    return total_inconsistencies, seq_f-seq_i
    """
    par_i = time.time()
    total_inconsistencies2 = inconsistency_finder.write_inconsistencies_parallel(
        inconsistency_fname=out_name+"para", num_workers=12, m_partitions=3)
    par_f = time.time()
    with open("total_inconsistency.txt", 'w') as f:
        f.write("Total number of Inconsistencies (Sequential) Found: " + str(total_inconsistencies) + "\n")
        f.write("Total number of Inconsistencies (Parallel) Found: " + str(total_inconsistencies2) + "\n")
        f.write("Sequential Time : " + str(seq_f-seq_i) + "\n")
        f.write("Parallel Time : " + str(par_f-par_i))

    print("Total number of Inconsistencies (Sequential) Found: " + str(total_inconsistencies))
    print("Total number of Inconsistencies (Parallel) Found: " + str(total_inconsistencies2))
    print("Sequential Time : " + str(seq_f-seq_i))
    print("Parallel Time : " + str(par_f-par_i))

    """
# ----------------------


def main():
    if len(sys.argv) != 3:
        print("USAGE: python revco.py <parsed blast file> <MEGAN_folder(MustHave=Species,Family,Class)>")
        exit(-1)
        
    print("############# REVCO STARTED #############")
    print("------- GREFCO STARTED -------")
    grefco = grefco_main()
    print("------- GREFCO FINISHED -------")
    print("------- TAXONOMY_LOAD STARTED -------")
    print("<-> Species")
    species_read_tax, species_contig_tax = taxonomy_load_main(tax_level="Species")
    print("<-> Family")
    family_read_tax, family_contig_tax = taxonomy_load_main(tax_level="Family")
    print("<-> Class")
    class_read_tax, class_contig_tax = taxonomy_load_main(tax_level="Class")
    print("------- TAXONOMY_LOAD FINISHED -------")
    print("------- INCONSISTENCY_FINDER STARTED -------")
    sp_n_inc, sp_time = inconsistency_finder_main(grefco, species_read_tax, species_contig_tax, "species_synth.inconsistency")
    fm_n_inc, fm_time = inconsistency_finder_main(grefco, family_read_tax, family_contig_tax, "family_synth.inconsistency")
    cl_n_inc, cl_time = inconsistency_finder_main(grefco, class_read_tax, class_contig_tax, "class_synth.inconsistency")
    with open("total_inconsistency.txt", 'w') as f:
        f.write("Taxon Level;Number of Inconsistencies;Time\n")
        f.write("Species;" + str(sp_n_inc) + ";" + str(sp_time) + "\n")
        f.write("Family;" + str(fm_n_inc) + ";" + str(fm_time) + "\n")
        f.write("Class;" + str(cl_n_inc) + ";" + str(cl_time) + "\n")
    print("------- INCONSISTENCY_FINDER FINISHED -------")
    print("############# REVCO FINISHED #############")


# ----------------------


if __name__ == "__main__":
    main()
