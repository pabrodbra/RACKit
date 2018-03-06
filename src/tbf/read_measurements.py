import sys

def main():
	if(len(sys.argv) != 2):
		print("***ERROR*** | USE: python read_measurement.py <fixed_reads_blast>"); quit()
	fix_blast_reads = sys.argv[1]
	tp=0
	tn=0
	fp=0
	fn=0
	total=0
	with open(fix_blast_reads, 'r') as f:
		for line in f:
			items = line.split(' ')
			r_id = items[0].rsplit('.',1)[0]
			#print("R_ID_ORIG :: " + r_id)
			r_match = items[1]
			#print("R_MATCH :: " + r_match)
			if (r_id==r_match and "RANDOMREAD" not in r_id):
				tp = tp + 1
			elif (r_id!=r_match and "RANDOMREAD" not in r_id):
				fp = fp + 1
			elif (r_id == r_match):
				fn = fn + 1
				#print(r_id + " :: " + r_match)
			else:
				tn = tn + 1
				#print(r_id + " :: " + r_match)
			total = total + 1
	spec = tn/float(tn+fp)
	sens = tp/float(tp+fn)
	correct = (tp + tn)/float(total)
	incorrect = (fp + fn)/float(total)
	print("*** READS ***")
	print("TP;FP;FN;TN;TotalRelationships;Specificity;Sensitivity;Correct;Incorrect")
	print(str(tp)+";"+str(fp)+";"+str(fn)+";"+str(tn)+";"+str(total)+";"+str(spec)+";"+str(sens)+";"+str(correct)+";"+str(incorrect))

#---------------

if __name__ == "__main__":
	main()
