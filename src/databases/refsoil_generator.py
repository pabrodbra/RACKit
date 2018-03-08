import sys
import requests
import csv

############################
# Obtain FASTA DB from NCBI ID from RefSoil_V1
############################

def main():
	fname = sys.argv[1]
	output_path = sys.argv[2]

	with open(fname, 'r') as csvf:
		reader = csv.reader(csvf, delimiter='\t')
		next(reader, None)

		counter = 0
		with open("ids_downloaded.log", 'w') as log_f:
			with open(output_path, 'w') as fo:
				for line in reader:
					ncbi_ids = line[2].split(',')

					for ncbi_id in ncbi_ids:

						ncbi_url = "http://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?tool=portal&sendto=on&log$=seqview&db=nuccore&dopt=fasta&val=%s&extrafeat=0&maxplex=1" % ncbi_id
						
						try:
							fasta_data = requests.get(ncbi_url, stream=True)
							for chunk in fasta_data.iter_content(chunk_size=1024*10):
								if chunk:
									text = chunk.decode('utf-8')
									fo.write(text)

							log_f.write(ncbi_id + "\n")

							counter = counter + 1
							if counter % 100 == 0:
								print("Genomes Downloaded: " + str(counter))
							
						except requests.exceptions.RequestException as e:
							print("ERROR | " + str(e))

						

# ---------------------

if __name__ == "__main__":
	if(len(sys.argv) != 3):
		print("ERROR | USE: python refsoil_generator.py <refsoil_txt> <output_DB>")
		quit()

	main()