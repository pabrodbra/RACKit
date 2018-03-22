from classes.grefco import GREFCO
from classes.taxaLoad import TaxaLoad
from classes.inconsistencyFinder import InconsistencyFinder
from classes.inconsistencySolver import InconsistencySolver
from classes.confusionMatrix import ConfusionMatrixCalculator
import sys
import time

DEBUG = True
OUTPUT_DIRECTORY = ""

def grefco_main(filtered_file_name, output_extension = ".grefco"):
    grefco = GREFCO(filtered_file_name)

    grefco.make_full_read_dictionary()
    grefco.make_full_contig_dictionary()

    if(DEBUG):
        GREFCO.output_dictionary(OUTPUT_DIRECTORY + "contig" + output_extension, grefco.contig_dictionary)
        GREFCO.output_dictionary(OUTPUT_DIRECTORY + "read" + output_extension, grefco.reads_dictionary)

    return grefco


def taxonomy_load_main(read_taxa_path, contigs_taxa_path, tax_level="Specie"):
    read_taxa_load = TaxaLoad(read_taxa_path, rank=tax_level)
    contig_taxa_load = TaxaLoad(contigs_taxa_path, rank=tax_level)

    if(DEBUG):
        read_taxa_load.save_taxonomy_dictionary(OUTPUT_DIRECTORY + tax_level + "read_category.csv")
        contig_taxa_load.save_taxonomy_dictionary(OUTPUT_DIRECTORY + tax_level + "contig_category.csv")
    
    return read_taxa_load, contig_taxa_load


def inconsistency_finder_main(grefco, read_tax, contig_tax, output = "inc_finder.out"):
    inconsistency_finder = InconsistencyFinder(grefco.reads_dictionary,
                                                       read_tax.get_taxonomy_dictionary(),
                                                       contig_tax.get_taxonomy_dictionary())

    seq_i = time.time()
    total_inconsistencies = inconsistency_finder.write_inconsistencies_sequential(inconsistency_fname=output)
    seq_f = time.time()
    total_time = seq_f-seq_i

    if(DEBUG):
        with open(OUTPUT_DIRECTORY + "total_inconsistency.txt", 'w') as f:
            f.write("Taxon Level;Number of Inconsistencies;Time\n")
            f.write("Species;" + str(total_inconsistencies) + ";" + str(total_time) + "\n")


def inconsistency_solver_main(read_taxa_path, contig_taxa_path, inconsistency_file, rank_filter, output = "inc_solver.out"):
    inc_solver = InconsistencySolver(read_taxa_path, contig_taxa_path, inconsistency_file)
    inc_solver.solve_inconsistencies(rank_filter, output)

def confusion_matrix_main(fixed_reads_blast, fixed_contigs_blast, grefco, grinder_ranks_path):
    confusion_matrix_calc = ConfusionMatrixCalculator(fixed_reads_blast, fixed_contigs_blast, grefco, grinder_ranks_path)

    tot_read = confusion_matrix_calc.read_confusion_matrix()
    tot_cont = confusion_matrix_calc.contig_confusion_matrix()
    confusion_matrix_calc.summarize_confusion_matrices(tot_read, tot_cont)
    confusion_matrix_calc.calculate_statistical_measurements()

    confusion_matrix_calc.write_confusion_matrix(OUTPUT_DIRECTORY + "rac-confusion_matrix.csv")
    confusion_matrix_calc.write_statistical_measurements(OUTPUT_DIRECTORY + "rac-statistical_measurements.csv")

# ----------------------

def main():
    if len(sys.argv) != 9:
        print("USAGE: python rackit.py <parsed_read_vs_contigs> <reads_taxon_path> <contigs_taxon_path> <rank_filter> <output_dir> <fixed_read_blast> <fixed_contig_blast> <grinder_ranks_path>")
        exit(-1)
        
    global OUTPUT_DIRECTORY
    global DEBUG

    DEBUG=True

    # Init Args
    filter_parsed_blast_arg = sys.argv[1]
    read_taxa_path_arg = sys.argv[2]
    contig_taxa_path_arg = sys.argv[3]
    rank_filter_arg = sys.argv[4]
    OUTPUT_DIRECTORY = sys.argv[5]
    fixed_reads_blast = sys.argv[6]
    fixed_contigs_blast = sys.argv[7]
    grinder_ranks_path = sys.argv[8]

    print("############# RACKIT.py STARTED #############")
    print("--- (1/5) GREFCO started...")
    grefco = grefco_main(filter_parsed_blast_arg)
    print("--! (1/5) GREFCO finished...")


    print("--- (2/5) TAXA_LOAD started...")
    print("<-> Species")
    species_read_tax, species_contig_tax = taxonomy_load_main(read_taxa_path_arg, contig_taxa_path_arg, tax_level="Specie")
    print("--- (2/5) TAXONOMY_LOAD finished...")


    print("--- (3/5) INCONSISTENCY_FINDER started...")
    inconsistency_finder_output = OUTPUT_DIRECTORY + "inc_finder.out"
    inconsistency_finder_main(grefco, species_read_tax, species_contig_tax, inconsistency_finder_output)
    print("--- (3/5) INCONSISTENCY_FINDER finished...")

    print("--- (4/5) INCONSISTENCY_SOLVER started...")
    inconsistency_solver_output = OUTPUT_DIRECTORY + "inc_solver.out"
    inconsistency_solver_main(read_taxa_path_arg, contig_taxa_path_arg, inconsistency_finder_output, rank_filter_arg, inconsistency_solver_output)
    print("--! (4/5) INCONSISTENCY_SOLVER finished...")

    print("--- (5/5) CONFUSION_MATRIX started...")
    confusion_matrix_main(fixed_reads_blast, fixed_contigs_blast, grefco, grinder_ranks_path)
    print("--! (5/5) CONFUSION_MATRIX finished...")

    print("############# REVCO FINISHED #############")

# ----------------------

if __name__ == "__main__":
    main()
