from classes.inconsistencySolver import InconsistencySolver

import sys
import time

# ----------------------
# python solve_inconsistency.py "C:\Users\Blinsky.Blinsk\Documents\- UMA\BitLab\Metagenome\Synthetic Dataset\SemiSynth\MEGAN\BLAST_ssd_reads-path.txt" "C:\Users\Blinsky.Blinsk\Documents\- UMA\BitLab\Metagenome\Synthetic Dataset\SemiSynth\MEGAN\BLAST_ssd_contigs-path.txt" "C:\Users\Blinsky.Blinsk\Documents\- UMA\BitLab\Metagenome\Synthetic Dataset\SemiSynth\results\species_synth.inconsistency" 10 ssd_solved

def main():
    if len(sys.argv) != 6:
    	print("ARGV length: " + str(len(sys.argv)))
    	print("USAGE: python revco.py <read_taxo_path> <contig_taxo_path> <inconsistencies> <taxo_rank_filter> <output>")
    	exit(-1)

    inc_solver = InconsistencySolver(sys.argv[1], sys.argv[2], sys.argv[3])
    inc_solver.solve_inconsistencies(sys.argv[4], sys.argv[5])


# ----------------------

if __name__ == "__main__":
    main()
