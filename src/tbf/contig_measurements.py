import sys
#from grefco import GREFCO

def contig_measurement(fix_blast_contigs, grefco_contig):
        tp=0
        tn=0
        fp=0
        fn=0
        total=0
        # LOAD CONTIG DICT
        grefco = GREFCO(grefco_contig)
        grefco.make_full_contig_dictionary()
        # CREATE REPRESENTATIVE DICTIONARY
        contig_most_representative = dict()
        for c_id in grefco.contig_dictionary.keys():
                temp_dict = {}
                asm_reads = grefco.contig_dictionary.get(c_id, None)
                for t_read in asm_reads:
                        r_id = t_read[0].rsplit('.',1)[0]
                        temp_dict[r_id] = temp_dict.get(r_id, 0) + 1
                        #if r_id not in temp_dict.keys():
                        #        temp_dict[r_id] = 1
                        #else:
                        #        temp_dict[r_id] = temp_dict[r_id] + 1
                max_v = max(temp_dict.values())
                tops = [k for k,v in temp_dict.items() if v == max_v]
                contig_most_representative[c_id] = tops;
        # GET MATRIX VALUES
        with open(fix_blast_contigs, 'r') as f:
                for line in f:
                        items = line.split(' ')
                        contig_id = items[0][1:]
                        contig_match = items[1][1:]
                        res = contig_most_representative.get(contig_id, None)
                        if res is not None:
                                for top_read in res:
                                        if (top_read == contig_match and "RANDOMREAD" not in contig_id):
                                                tp = tp + 1
                                        elif (top_read != contig_match and "RANDOMREAD" not in contig_id):
                                                fp = fp + 1
                                        elif (top_read == contig_match):
                                                fn = fn + 1
                                        	#print(r_id + " :: " + r_match)
                                        else:
                                                tn = tn + 1
                                        #print(r_id + " :: " + r_match)
                                        total = total + 1
        # CALCULATE
        spec = tn/float(tn+fp)
        sens = tp/float(tp+fn)
        correct = (tp + tn)/float(total)
        incorrect = (fp + fn)/float(total)
        # PRINT
        print("TP;FP;FN;TN;TotalRelationships;Specificity;Sensitivity;Correct;Incorrect")
        print(str(tp)+";"+str(fp)+";"+str(fn)+";"+str(tn)+";"+str(total)+";"+str(spec)+";"+str(sens)+";"+str(correct)+";"+str(incorrect))

def main():
        if(len(sys.argv) != 3):
                print("***ERROR*** | USE: python read_measurement.py <fixed_contigs_blast> <contig_read_filter_blast>"); quit()
        fix_blast_contigs = sys.argv[1]
        contig_grefco = sys.argv[2]
        contig_measurement(fix_blast_contigs, contig_grefco)

if __name__ == '__main__':
        main()
