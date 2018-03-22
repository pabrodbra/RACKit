from classes.grefco import GREFCO

def make_confusion_matrix_dict():
    tp = tn = fp = fn = 0
    confusion_matrix = {"tp": tp, "tn": tn, "fp": fp, "fn": fn}
    return confusion_matrix

def make_statistical_measurements_dict(cm):
    tp = cm["tp"]; fp = cm["fp"]; tn = cm["tn"]; fn = cm["fn"]
    sens = tn/float(tn+fp)
    spec = tp/float(tp+fn)
    fallout = fp/float(fp+tn)
    precision = tp/float(tp+fp)
    acc = (tp + tn)/float(tp + fp + fn + tn)
    measures = {"sensitivity": sens, "specificity": spec, "fallout": fallout, "precision": precision, "accuracy": acc}
    return measures

def confusion_matrix_to_string(cm):
    ret = str(cm["tp"]) + ";" + str(cm["fp"]) + ";" + str(cm["tn"]) + ";" + str(cm["fn"])
    return ret

def statistical_measurements_to_string(sm):
    ret = str(sm["sensitivity"]) + ";" + str(sm["specificity"]) + ";" + str(sm["fallout"]) + ";" + str(sm["precision"]) + ";" + str(sm["accuracy"])
    return ret

class ConfusionMatrixCalculator(object):
    def __init__(self, fixed_read_blast, fixed_contig_blast, grefco, grinder_rank):
        self.fixed_reads = fixed_read_blast
        self.fixed_contigs = fixed_contig_blast
        self.grefco = grefco
        self.reads_all_confusion_matrix = {}
        self.contigs_all_confusion_matrix = {}

        with open(grinder_rank, 'r') as f:
            header = f.readline()
            for line in f:
                items = line.split('\t')
                self.reads_all_confusion_matrix[items[1]] = make_confusion_matrix_dict()
                self.contigs_all_confusion_matrix[items[1]] = make_confusion_matrix_dict()

    def read_confusion_matrix(self):
        total = 0
        with open(self.fixed_reads, 'r') as f:
            # >5-NC_014963.1 >NC_014963.1.1  5095226 1;0;90;90;100;0;0;+-;1;90;3586622;3586533
            for line in f:
                items = line.split(' ')
                r_id = items[0].split('-')[1] # >5-NC_014963.1 -> NC_014963.1
                r_match = items[1].rsplit('.', 1)[0][1:] # >NC_014963.1.1 -> NC_014963.1
                
                for specie, s_cm in self.reads_all_confusion_matrix.items():
                    if (r_match == specie):
                        if(r_match == r_id):
                           s_cm["tp"] += 1
                        else:
                            s_cm["fp"] += 1
                    
                    if (r_match != specie):
                        if(r_match == r_id):
                            s_cm["tn"] += 1
                        else:
                            s_cm["fn"] += 1
                    total += 1

        return total

    def obtain_contigs_representative_reads(self):
        # Make representative dictionary of which reads make most of a contigs
        contig_most_representative = {}
        for c_id in self.grefco.contig_dictionary.keys():
                temp_dict = {}
                asm_reads = self.grefco.contig_dictionary.get(c_id, None)
                for t_read in asm_reads:
                        r_id = t_read[0].rsplit('-',1)[1]
                        temp_dict[r_id] = temp_dict.get(r_id, 0) + 1
                max_v = max(temp_dict.values())
                tops = [k for k,v in temp_dict.items() if v == max_v]
                contig_most_representative[c_id] = tops

        return contig_most_representative


    def contig_confusion_matrix(self):
        contig_most_representative = self.obtain_contigs_representative_reads()
        total = 0
        with open(self.fixed_contigs, 'r') as f:
            # >k67_1 flag=1 multi=2.0000 len=207 >NC_014963.1.1  5095226 1;0;207;207;100;0;0;+-;1;207;4289335;4289129
            for line in f:
                items = line.split(' ')
                contig_id = items[0][1:] # >k67_1 -> k67_1
                contig_match = items[4][1:].rsplit('.',1)[0] # >NC_014963.1.1 -> NC_014963.1
                # print(contig_id); print(contig_match)
                representatives = contig_most_representative.get(contig_id, None) #; print(representatives); input()
                if representatives is not None:
                    for top_read in representatives:
                        for specie, s_cm in self.contigs_all_confusion_matrix.items():
                            if (contig_match == specie):
                                if(contig_match == top_read):
                                    s_cm["tp"] += 1
                                else:
                                    s_cm["fp"] += 1
                            
                            if (contig_match != specie):
                                if(contig_match == top_read):
                                    s_cm["tn"] += 1
                                else:
                                    s_cm["fn"] += 1
                            total += 1

        return total

    def summarize_confusion_matrices(self, total_reads, total_contigs):
        reads_final_cm = make_confusion_matrix_dict(); print(self.reads_all_confusion_matrix)
        contigs_final_cm = make_confusion_matrix_dict(); print(self.contigs_all_confusion_matrix)

        for specie, confusion_matrix in self.reads_all_confusion_matrix.items():
            reads_final_cm["tp"] += confusion_matrix["tp"]/float(total_reads)
            reads_final_cm["tn"] += confusion_matrix["tn"]/float(total_reads)
            reads_final_cm["fp"] += confusion_matrix["fp"]/float(total_reads)
            reads_final_cm["fn"] += confusion_matrix["fn"]/float(total_reads)

        for specie, confusion_matrix in self.contigs_all_confusion_matrix.items():
            contigs_final_cm["tp"] += confusion_matrix["tp"]/float(total_contigs)
            contigs_final_cm["tn"] += confusion_matrix["tn"]/float(total_contigs)
            contigs_final_cm["fp"] += confusion_matrix["fp"]/float(total_contigs)
            contigs_final_cm["fn"] += confusion_matrix["fn"]/float(total_contigs)

        self.final_read_cm = reads_final_cm; print(reads_final_cm)
        self.final_contig_cm = contigs_final_cm

    def calculate_statistical_measurements(self):
        read_sm = make_statistical_measurements_dict(self.final_read_cm)
        contig_sm = make_statistical_measurements_dict(self.final_contig_cm)

        self.final_read_sm = read_sm
        self.final_contig_sm = contig_sm

    def write_confusion_matrix(self, output="rac-confusion_matrix.csv"):
        with open(output, 'w') as f:
            f.write("TruePositive;FalsePositive;TrueNegative;FalseNegative")
            f.write(confusion_matrix_to_string(self.final_read_cm))
            f.write(confusion_matrix_to_string(self.final_contig_cm))

    def write_statistical_measurements(self, output="rac-statistical_measurements.csv"):
        with open(output, 'w') as f:
            f.write("Sensitivity;Specificity;Fallout;Precision;Accuracy")
            f.write(statistical_measurements_to_string(self.final_read_sm))
            f.write(statistical_measurements_to_string(self.final_contig_sm))