from classes.grefco import GREFCO
from classes.taxaLoad import TaxaLoad
from classes.inconsistencyFinder import InconsistencyFinder

import sys
import time

DEBUG = True
OUTPUT_DIRECTORY = ""

def grefco_main(filtered_file_name, output = ".grefco"):
    grefco = GREFCO(filtered_file_name)

    grefco.make_full_read_dictionary()
    grefco.make_full_contig_dictionary()

    if(DEBUG):
        GREFCO.output_dictionary(OUTPUT_DIRECTORY + "contig" + output, grefco.contig_dictionary)
        GREFCO.output_dictionary(OUTPUT_DIRECTORY + "read" + output, grefco.reads_dictionary)

    return grefco


def taxonomy_load_main(read_taxa_path, contigs_taxa_path, tax_level="Species"):
    read_taxa_load = TaxaLoad(read_taxa_path, rank=tax_level)
    contig_taxa_load = TaxaLoad(contigs_taxa_path, rank=tax_level)

    if(DEBUG):
        read_taxa_load.save_taxonomy_dictionary("species-read_category.csv")
        contig_taxa_load.save_taxonomy_dictionary("species-contig_category.csv")
    
    return read_taxa_load, contig_taxa_load


def inconsistency_finder_main(grefco, read_tax, contig_tax, output = "inc_finder.out"):
    inconsistency_finder = InconsistencyFinder(grefco.reads_dictionary,
                                                       read_tax.get_taxonomy_dictionary(),
                                                       contig_tax.get_taxonomy_dictionary())

    seq_i = time.time()
    total_inconsistencies = inconsistency_finder.write_inconsistencies_sequential(
        inconsistency_fname=output)
    seq_f = time.time()
    total_time = seq_f-seq_i

    print("Total number of Inconsistencies Found: " + str(total_inconsistencies))
    print("Sequential Time: " + str(total_time))

    if(DEBUG):
        with open("total_inconsistency.txt", 'w') as f:
            f.write("Taxon Level;Number of Inconsistencies;Time\n")
            f.write("Species;" + str(total_inconsistencies) + ";" + str(total_time) + "\n")


def inconsistency_solver_main(output = "inc_solver.out"):
    return 0

def statistical_measurements_main(output = "statistic.out"):
    return 0

# ----------------------

def main():
    if len(sys.argv) != 3:
        print("USAGE: python rackit.py <parsed_read_vs_contigs> <reads_taxon_path> <contigs_taxon_path> <rank_filter> <output_dir> (inconsitency)]")
        exit(-1)
        
    global OUTPUT_DIRECTORY
    OUTPUT_DIRECTORY = sys.argv[5]


    print("############# RACKIT.py STARTED #############")
    print("--- (1/5) GREFCO started...")
    grefco = grefco_main(sys.argv[1])
    print("--! (1/5) GREFCO finished...")


    print("--- (2/5) TAXA_LOAD started...")
    print("<-> Species")
    species_read_tax, species_contig_tax = taxonomy_load_main(sys.argv[2], sys.argv[3], tax_level="Species")
    print("--- (2/5) TAXONOMY_LOAD finished...")


    print("--- (3/5) INCONSISTENCY_FINDER started...")
    inconsistency_finder_main(grefco, species_read_tax, species_contig_tax, OUTPUT_DIRECTORY + "found_inconsistencies.inc")
    print("--- (3/5) INCONSISTENCY_FINDER finished...")

    print("--- (4/5) INCONSISTENCY_SOLVER started...")
    
    print("--! (4/5) INCONSISTENCY_SOLVER finished...")


    print("--- (5/5) STATISTICAL_MEASUREMENTS started...")

    print("--! (5/5) STATISTICAL_MEASUREMENTS finished...")

    print("############# REVCO FINISHED #############")

# ----------------------

if __name__ == "__main__":
    main()
