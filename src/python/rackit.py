from classes.grefco import GREFCO
from classes.taxonomyLoad import TaxonomyLoad
from classes.inconsistencyFinder import InconsistencyFinder

import sys
import time

def grefco_main():
    return 0

def taxonomy_load_main(tax_level):
    return 0,0

def inconsistency_finder_main(grefco, read_tax, contig_tax, output_name):
    return 0, 0

def inconsistency_solver_main():
    return 0

def statistical_measurements_main():
    return 0

# ----------------------

def main():
    if len(sys.argv) != 3:
        print("USAGE: python rackit.py  [<parsed_read_vs_contigs> <Taxa_folder(Reads,Contigs inside)>] [<reads_taxon_path> <contigs_taxon_path> <rank_filter> <output> (inconsitency)]")
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
