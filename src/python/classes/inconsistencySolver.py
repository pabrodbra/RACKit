from __future__ import print_function, division

class InconsistencySolver(object):
    def __init__(self, reads_taxo_path, contigs_taxo_path, inconsistencies_file):
        self.inconsistency_file = inconsistencies_file
        self.reads_path_dict = dict()
        self.contigs_path_dict = dict()

        with open(reads_taxo_path, 'r') as f:
            for line in f:
                items = line.split(",")
                path = items[1].replace('"','').split(";")[:-1]
                self.reads_path_dict[items[0]] = path

        with open(contigs_taxo_path, 'r') as f:
            for line in f:
                items = line.split(",")
                path = items[1].replace('"','').split(";")[:-1]
                self.contigs_path_dict[items[0]] = path

    def solve_inconsistencies(self, rank_filter, output_file):
        weak_counter_reads = 0
        weak_counter_contigs = 0
        hard_counter = 0
        line_counter = 0
        filter_rank = int(rank_filter)

        with open(output_file, 'w') as fout:
            fout.write("InconsistencyLine,ReadID,ReadTopRank,ReadTopTaxo,ContigID,ContigTopRank,ContigTopTaxo,Solvable,RankSolved,Solution,Type\n")
            with open(self.inconsistency_file, 'r') as f:
                for line in f:
                    line_counter += 1
                    items = line.split(",")
                    read_id = items[0]
                    contig_id = items[2]
                    read_top_rank = contig_top_rank = -1
                    read_top_taxo = contig_top_taxo = solvable = type_inconsistency = ""

                    # Get Top Rank + Taxo
                    if read_id in self.reads_path_dict:
                        read_top_rank = len(self.reads_path_dict[read_id])
                        read_top_taxo = self.reads_path_dict[read_id][read_top_rank-1]
                    if contig_id in self.contigs_path_dict:
                        contig_top_rank = len(self.contigs_path_dict[contig_id])
                        contig_top_taxo = self.contigs_path_dict[contig_id][contig_top_rank-1]

                    # Get Solvable + RanksToSolve + RankSolved + Type
                    if read_top_rank != -1 and contig_top_rank != -1:   # Can be solved
                        solvable = "Y"

                        read_path = self.reads_path_dict[read_id]
                        contig_path = self.contigs_path_dict[contig_id]

                        solved_ranks = list(set(read_path).intersection(contig_path))
                        ranks_to_solve = len(solved_ranks)
                        solution = read_path[ranks_to_solve-1]

                        # Get Type
                        if read_top_rank == contig_top_rank:
                            hard_counter += 1
                            type_inconsistency = "H"
                        else:
                            if read_top_rank > contig_top_rank:
                                weak_counter_contigs += 1
                                type_inconsistency = "WC"
                            else:
                                weak_counter_reads += 1
                                type_inconsistency = "WR"

                    else:   # Cannot be solved
                        solvable = "N"

                        ranks_to_solve = -1
                        solution = "None"

                        # Get Type
                        if read_top_rank == -1:
                            type_inconsistency = "WR"
                        else:
                            type_inconsistency = "WC"

                    # Write output
                    if ranks_to_solve > -1 and ranks_to_solve < filter_rank:
                        fout.write(str(line_counter)
                                   + "," +  read_id + "," +  str(read_top_rank) + "," +  read_top_taxo
                                   + "," +  contig_id + "," +  str(contig_top_rank) + "," +  contig_top_taxo
                                   + "," +  solvable + "," + str(ranks_to_solve) + "," +  solution
                                   + "," + type_inconsistency + "\n")
