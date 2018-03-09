from __future__ import print_function
import sys
import subprocess, shlex

# Print function for stderr and quit
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)
	quit()

############################
# Obtain FASTA from NCBI ID
############################

def main():
	source_fname = sys.argv[1]
	base_metagenome = sys.argv[2]
	read_length = sys.argv[3]
	output_path = sys.argv[4]

	with open(source_fname, 'r') as f:
		for line in f:
			items = line.replace('\n','').split(',')
			n_reads = items[0]
			name_org = items[1]
			pa_out = output_path + name_org.split('/')[-1].split('.')[0] + ".pa"
			fasta_out = output_path + "reads_" + name_org.split('/')[-1].split('.')[0] + ".fasta"
			fastq_out = output_path + "reads_" + name_org.split('/')[-1].split('.')[0] + ".fastq"

			pa_args = name_org + " high"
			gr_args = "-r " + name_org + " -a " + pa_out + " -o " + fasta_out + " -t " + n_reads + " -l " + read_length
			te_args = "-i " + base_metagenome + " -f " + fasta_out + " -o " + fastq_out + " -r 0 -q 0 -m 0"

			f_pa = open(pa_out, 'w')

			pa = shlex.split("/home/estebanpw/software/BEAR/BEAR-master/scripts/parametric_abundance.pl " + pa_args)
			gr = shlex.split("/home/estebanpw/software/BEAR/BEAR-master/scripts/generate_reads.py " + gr_args)

			p = subprocess.Popen(pa, stdout=f_pa)
			p.wait()
			f_pa.close()
			p = subprocess.Popen(gr)
			p.wait()
			'''
			te = shlex.split("/home/pablorod/data/metagenomics/BEAR/trim_edit.pl " + te_args)
			fq2fa = shlex.split("/home/pablorod/software/scripts/fastq2fasta.sh " + fastq_out)
			p = subprocess.Popen(te)
			p.wait()
			p = subprocess.Popen(fq2fa)
			p.wait()
			'''

			print('READS of %s generated succesfully!' % name_org)


# ---------------------

if __name__ == "__main__":
	if(len(sys.argv) != 5):
		# python generate_reads.py /home/pablorod/data/metagenomics/semisynth/ /home/pablorod/data/metagenomics/semisynth/mycoplasma/mycoplasma_7422_FKHDS0V01.fasta 250 ./
		eprint("ERROR | USE: python generate_reads.py <source_fasta.csv(read_quantity,fasta)> <base_metagenome> <r_length> <output_folder>")

	main()

'''
BASE_METAGENOME=mycoplasma_7422_FKHDS0V01.fasta

@BEAR FOLDER
/home/estebanpw/software/BEAR/BEAR-master/scripts/parametric_abundance.pl [1] high > [1]+".pa"
/home/estebanpw/software/BEAR/BEAR-master/scripts/generate_reads.py -r [1] -a [1]+".pa" -o "reads_"+[1]+".fasta" -t [0] -l 250
/home/pablorod/data/metagenomics/BEAR/trim_edit.pl -i BASE_METAGENOME -f "reads_"+[1]+".fasta" -o "reads_"+[1]+".fastq" -r 0 -q 0 -m 0
fastq2fasta "reads_"+[1]+".fastq
'''