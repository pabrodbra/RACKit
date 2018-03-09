import sys
import math
#from PIL import Image, ImageFont, ImageDraw

def load_seq_file(fname):
	sid = ""
	seq = ""
	with open(fname,'r') as f:
		content = f.readlines()
		sid = content[0]
		seq = "".join(content[1:]).replace('\n','')

		return sid, seq, len(seq)

def RAC_coverage_comparison(s_len, read_seq, contig_seq, max_hits):
	read_only = 0
	contig_only = 0
	both_yes = 0
	both_no = 0

	read_hits = 0
	contig_hits = 0

	read_mat = list()
	contig_mat = list()

	for i in range(0, s_len-1, 1):
		if read_seq[i]==contig_seq[i]:
			if read_seq[i] == '1':
				both_yes+=1
			else:
				both_no+=1
		else:
			if read_seq[i] == '1':
				read_only+=1
				read_hits+=1
			else:
				contig_only+=1
				contig_hits+=1

		if (i+1) % max_hits == 0:
			read_mat.append(read_hits)
			contig_mat.append(contig_hits)
			read_hits = contig_hits = 0


	return read_only, contig_only, both_yes, both_no, read_mat, contig_mat

"""
def matrix_to_png(read_mat, contig_mat, max_hits, img_dim):
	MAX_COLOR = 255
	FACTOR = MAX_COLOR / max_hits

	img = Image.new("RGB", (img_dim, img_dim+100))

	draw = ImageDraw.Draw(img)
	font_large = ImageFont.load_default()

	for i in range(0, len(read_mat)-1, 1):
		curr_y = math.floor(i / img_dim)
		curr_x = i % img_dim
		r_val = math.floor(read_mat[i]*FACTOR)
		c_val = math.floor(contig_mat[i]*FACTOR)

		img.putpixel((curr_x, curr_y), (r_val, 0, c_val))

	draw.rectangle(((img_dim/10 - 10, img_dim+40), (img_dim/10 + 10, img_dim+60)), fill="red")
	draw.rectangle(((3*img_dim/10 - 10, img_dim+40), (3*img_dim/10 + 10, img_dim+60)), fill="blue")
	draw.text((img_dim/10, img_dim+50), "Reads", (255,255,255), font=font_large)
	draw.text((3*img_dim/10, img_dim+50), "Contig", (255,255,255), font=font_large)
	draw.text((5*img_dim/10 + 50, img_dim+50), "Both", (255,255,255), font=font_large)

	img.save("RAC_coverage.png")
"""

def save_matrix(matrix, fout, max_hits):
	with open(fout, 'w') as f:
		f.write("Pixel;Hits<" + max_hits)
		for i in range(0,len(matrix)-1,1):
			f.write(str(i) + ";" + str(matrix[i]))


############################
# Create Scaled Matrix from 2 '.seq' files (from uniseqDBCoverage) + Draw + Coverage Comparison
############################

def main():
	read_seq_fname = sys.argv[1]
	contig_seq_fname = sys.argv[2]
	dim = 1024

	r_id, r_seq, r_len = load_seq_file(read_seq_fname)
	c_id, c_seq, c_len = load_seq_file(contig_seq_fname)

	max_hits_per_pixel = math.ceil(r_len / (dim * dim))

	print(max_hits_per_pixel)

	r_only, c_only, b_yes, b_no, r_mat, c_mat = RAC_coverage_comparison(r_len, r_seq, c_seq, max_hits_per_pixel)


	print("ReadOnly;ContigOnly;BothYes;BothNo")
	print(str(r_only) + ";" + str(c_only) + ";" + str(b_yes) + ";" + str(b_no))
	print("SUM: " + str(r_only + c_only + b_yes + b_no))

	with open("coverage_comparison.cc", 'w') as f:
		f.write("ReadOnly;ContigOnly;BothYes;BothNo\n")
		f.write(str(r_only) + ";" + str(c_only) + ";" + str(b_yes) + ";" + str(b_no))
		
	#matrix_to_png(r_mat, c_mat, max_hits_per_pixel, dim)
	save_matrix(r_mat, "Reads"+str(dim)+".csv")
	save_matrix(r_mat, "Contigs"+str(dim)+".csv")

# ---------------------

if __name__ == "__main__":
	if(len(sys.argv) != 3):
		print("ERROR | USE: python matrixSeqCompare.py <1.seq> <2.seq>"); quit()

	main()