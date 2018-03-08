from __future__ import print_function
import sys
import requests

# Print function for stderr and quit
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)
	quit()

############################
# Obtain FASTA from NCBI ID
############################

def main():
	fname = sys.argv[1]
	output_path = sys.argv[2]

	with open(fname, 'r') as f:
		for line in f:
			ncbi_id = line.replace('\n','').split(',')
			ncbi_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&val=%s&extrafeat=0&maxplex=1" % ncbi_id[0]
			
			try:
				fasta_data = requests.get(ncbi_url, stream=True)
			except requests.exceptions.RequestException as e:
				eprint("ERROR | " + str(e))

			f = open(output_path + ncbi_id[1], 'w')
			for chunk in fasta_data.iter_content(chunk_size=1024*10):
				if chunk:
					text = chunk.decode('utf-8')
					f.write(text)
			
			f.close()
			print('FASTA of %s downloaded succesfully!' % ncbi_id[0])


# ---------------------

if __name__ == "__main__":
	if(len(sys.argv) != 3):
		eprint("ERROR | USE: python nuccore_fasta.py <NCBI_IDs_file(id,fname)> <output_folder>")

	main()